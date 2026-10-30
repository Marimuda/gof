#!/usr/bin/env python3
"""Run compiler-backed checks on a frozen package using its actual install log.

Requires Python 3 and matching Clang tooling (defaults to version 19). No
dependencies are downloaded, source is not rewritten, and fixes are not applied.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", required=True, type=Path, help="frozen package source")
parser.add_argument("--install-log", required=True, type=Path, help="Clang R CMD INSTALL log")
parser.add_argument("--output", required=True, type=Path, help="new evidence directory")
parser.add_argument("--jobs", type=int, default=2)
parser.add_argument("--clang-tidy", default="clang-tidy-19")
parser.add_argument("--clang-check", default="clang-check-19")
parser.add_argument("--clang-format", default="clang-format-19")
args = parser.parse_args()
if args.jobs < 1:
    parser.error("--jobs must be positive")
source = args.source.resolve(strict=True)
install_log = args.install_log.resolve(strict=True)
output = args.output.resolve()
if output.exists() or output == source or source in output.parents:
    parser.error("choose a new output directory outside the package source")
database = []
for line in install_log.read_text().splitlines():
    if "clang" not in line or " -c " not in line:
        continue
    command = shlex.split(line)
    if "-c" not in command or "clang" not in Path(command[0]).name:
        continue
    unit = Path(command[command.index("-c") + 1]).name
    if not unit.endswith(".cpp") or not (source / "src" / unit).is_file():
        continue
    # R installs archives through a temporary source directory. Rebase the
    # command to the retained, byte-identical snapshot, preserving every flag
    # and relative include path. The compile database records this mapping.
    command[command.index("-c") + 1] = str(source / "src" / unit)
    database.append({"directory": str(source / "src"), "file": str(source / "src" / unit),
                     "arguments": command})
expected = {str(p) for p in (source / "src").glob("*.cpp")}
if {entry["file"] for entry in database} != expected:
    parser.error("install log does not cover every frozen translation unit; wait for the "
                 "Clang installation to finish and use its complete compiler log")
output.mkdir(parents=True)
(output / "logs").mkdir()
(output / "compile_commands.json").write_text(json.dumps(database, indent=2) + "\n")
owned = [Path(entry["file"]) for entry in database if Path(entry["file"]).name != "RcppExports.cpp"]
formatted = sorted(owned + list((source / "inst" / "include").rglob("*.hpp"))
                   + list((source / "tests" / "native").glob("*.cpp")))
inputs = sorted((source / "src").glob("*.cpp")) + sorted(
    p for p in (source / "inst" / "include").rglob("*") if p.suffix in (".h", ".hpp"))
inputs += [source / ".clang-format", source / ".clang-tidy"]
inputs += sorted((source / "tests" / "native").glob("*.cpp"))
hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
diagnostic = re.compile(r"^(.+?):(\d+):(\d+): (warning|error): (.*?)(?: \[([^]]+)\])?$")
records = []

def run(label, command, cwd=None):
    logfile = output / "logs" / (label + ".log")
    with logfile.open("w") as handle:
        result = subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=False)
    findings = []
    for line in logfile.read_text().splitlines():
        match = diagnostic.match(line)
        if match:
            findings.append(dict(zip(("file", "line", "column", "level", "message", "check"), match.groups())))
    return {"label": label, "command": command, "directory": str(cwd) if cwd else None,
            "exit_code": result.returncode,
            "log": str(logfile), "findings": findings}

for name, executable in (("tidy", args.clang_tidy), ("check", args.clang_check), ("format", args.clang_format)):
    records.append(run("version-" + name, [executable, "--version"]))
records.append(run("format", [args.clang_format, "--dry-run", "--Werror", *map(str, formatted)]))

# Check each owned header in isolation using the actual R-backed compiler flags.
header_command = list(database[0]["arguments"])
for option in ("-c", "-o"):
    index = header_command.index(option)
    del header_command[index:index + 2]
for header in sorted((source / "inst" / "include").rglob("*.hpp")):
    relative = header.relative_to(source / "inst" / "include")
    label = "header-" + str(relative).replace("/", "-").replace(".", "-")
    probe = output / (label + ".cpp")
    probe.write_text("#include <" + str(relative) + ">\n")
    records.append(run(label, header_command + ["-fsyntax-only", str(probe)], cwd=source / "src"))

def analyze(unit):
    name = unit.stem
    tidy = run("tidy-" + name, [args.clang_tidy, "-p", str(output),
        "--config-file=" + str(source / ".clang-tidy"),
        "--export-fixes=" + str(output / ("tidy-" + name + ".yaml")), str(unit)])
    analyzer = run("analyzer-" + name, [args.clang_check, "-p", str(output), "--analyze",
        "--extra-arg=-Xclang", "--extra-arg=-analyzer-output=plist-multi-file",
        "--analyzer-output-path=" + str(output / (name + ".plist")), str(unit)])
    return [tidy, analyzer]

with ThreadPoolExecutor(max_workers=args.jobs) as pool:
    for pair in pool.map(analyze, owned):
        records.extend(pair)
unchanged = all(hashlib.sha256((source / p).read_bytes()).hexdigest() == h for p, h in hashes.items())
failed = [r["label"] for r in records if r["exit_code"] != 0]
blocking_findings = [dict(record=r["label"], **finding)
    for r in records if r["label"].startswith("analyzer-")
    for finding in r["findings"]]
# clang-check may exit zero after reporting an analyzer warning. Treat its
# actual diagnostics as failures as well; a process return code is insufficient.
manifest = {"source": str(source), "install_log": str(install_log),
    "install_log_sha256": hashlib.sha256(install_log.read_bytes()).hexdigest(),
    "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "source_unchanged": unchanged,
    "compile_database_units": len(database), "owned_units_analyzed": len(owned),
    "generated_binding_exclusion": "RcppExports.cpp; generated headers excluded by owned-header filter",
    "input_sha256": hashes, "records": records, "failed_commands": failed,
    "blocking_analyzer_findings": blocking_findings,
    "note": "Review all findings. Advisory warnings are retained; exit status alone is not a clean-analysis claim."}
(output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"units": len(owned), "failed_commands": failed, "source_unchanged": unchanged,
                  "blocking_analyzer_findings": len(blocking_findings),
                  "displayed_findings": sum(len(r["findings"]) for r in records)}, indent=2))
if failed or blocking_findings or not unchanged:
    raise SystemExit(1)
