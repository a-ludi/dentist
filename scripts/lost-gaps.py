import argparse
import json
import re
import sys
from contextlib import contextmanager
from pathlib import Path

OnlyFlag = ("spanning", "extending", "both")
JoinPolicy = ("scaffoldGaps", "scaffolds", "contigs")


def get_error_cat(event):
    error = event.error
    if "consensus" in error:
        return "consensus failed"
    elif "no valid reference read found" in error:
        return "no valid reference read found"
    elif "empty pileup alignment" in error:
        return "reads in pile do not properly align to each other"
    elif error == "could not find a common trace point":
        return "read cropping failed"
    else:
        return "other"


def load_multi_json(fp, **kwargs):
    decoder = json.JSONDecoder(**kwargs)

    for i, line in enumerate(fp, 1):
        try:
            yield decoder.decode(line)
        except json.JSONDecodeError as e:
            raise ValueError(f"{e.msg}: {fp.name}, line {i}, column {e.pos}")


def total_len(mapping, *keys):
    tot = 0
    for key in keys:
        tot += len(mapping.get(key, []))
    return tot


def group_by(items, key):
    groups = namespace()

    for item in items:
        _key = key(item)
        _value = groups.setdefault(_key, list())
        _value.append(item)

    return groups


class namespace(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __getattr__(self, key):
        return self[key]


class Reporter:
    def __init__(self, log_dir, *, verbose):
        self.log_dir = Path(log_dir)
        self.verbose = verbose
        self._indent = 0

    def __call__(self, file=sys.stdout):
        self.collect_opts_events()
        self.events_by_reason = group_by(self.events, lambda e: e.reason)

        self.report_collect_phase()
        self.report_process_phase()
        self.report_output_phase()
        self.report_unhandled_events()

    @contextmanager
    def open_logs(self):
        logs = list()
        try:
            for pattern in ("collect.log", "process.*.log", "*-output.log"):
                for log_path in self.log_dir.glob(pattern):
                    logs.append(log_path.open())
            yield logs
        finally:
            for log in logs:
                log.close()
            del logs

    def collect_opts_events(self):
        raw_opts = namespace()
        events = list()

        with self.open_logs() as logs:
            for log in logs:
                for log_entry in load_multi_json(log, object_hook=namespace):
                    if "executableVersion" in log_entry:
                        raw_opts.update(log_entry)
                    elif "event" in log_entry and log_entry["event"] in [
                        "pileUpSkipped",
                        "insertionSkipped",
                    ]:
                        log_entry.setdefault("reason", "only")
                        events.append(log_entry)

        self.events = events
        self.opts = namespace(
            only=OnlyFlag[raw_opts.get("onlyFlag", 1) - 1],
            joinPolicy=JoinPolicy[raw_opts.get("joinPolicy", 0)],
            minSpanningReads=raw_opts.get("minSpanningReads", 3),
            minExtensionLength=raw_opts.get("minExtensionLength", 100),
            maxInsertionError=raw_opts.get("maxInsertionError", 0.1),
        )

    def cli_opt(self, opt, quote=True):
        if quote:
            return f"`{self.cli_opt(opt, quote=False)}`"
        cli_name = re.sub(r"([a-z])([A-Z])", lambda m: f"{m[1]}-{m[2].lower()}", opt)
        return f"--{cli_name}={self.opts[opt]}"

    def gap(self, thing):
        _gap = "-".join(str(cid) for cid in thing.contigIds)
        if "length" in thing:
            _gap += f" ({thing.length} reads)"
        return _gap

    def li(self, item):
        indent = self._indent * "    "
        print(f"{indent}- {item}")

    @contextmanager
    def indent(self):
        self._indent += 1
        yield
        self._indent -= 1

    def report_collect_phase(self):
        tot_lost = total_len(
            self.events_by_reason, "minSpanningReads", "scaffoldingConflict"
        )
        self.li(f"lost {tot_lost} in `collect` phase")

        with self.indent():
            self.report_min_spanning_reads()
            self.report_scaffolding_conflict()

    def report_min_spanning_reads(self):
        events = self.events_by_reason.pop("minSpanningReads", [])
        self.li(
            f"lost {len(events)} gap(s) because of insufficient number of "
            f"spanning reads ({self.cli_opt('minSpanningReads')})",
        )
        with self.indent():
            for event in events:
                self.li(f"skipped {self.gap(event.pileUp)}")

    def report_scaffolding_conflict(self):
        events = self.events_by_reason.pop("scaffoldingConflict", [])
        self.li(
            f"lost {len(events)} gap(s) because a scaffolding conflict was detected"
        )
        with self.indent():
            for event in events:
                gaps = (self.gap(pile_up) for pile_up in event.pileUps)
                self.li(f"conflicting gap closings: {', '.join(gaps)}")

    def report_process_phase(self):
        tot_lost = total_len(self.events_by_reason, "only", "error")
        self.li(f"lost {tot_lost} in `process` phase")

        with self.indent():
            self.report_error()
            self.report_only()

    def report_error(self):
        events = self.events_by_reason.pop("error", [])
        self.li(f"skipped {len(events)} read pile ups because of errors")
        by_error_cat = group_by(events, get_error_cat)
        with self.indent():
            for error_cat, event_list in by_error_cat.items():
                self.li(f"{error_cat} ({len(event_list)} times)")
                if self.verbose:
                    with self.indent():
                        for event in event_list:
                            self.li(f"{self.gap(event.pileUp)} failed: {event.error}")

    def report_only(self):
        events = self.events_by_reason.pop("only", [])
        self.li(
            f"skipped {len(events)} read pile ups because of {self.cli_opt('only')}",
        )
        if self.verbose:
            by_pile_up_type = group_by(events, lambda e: e.pileUp.type)
            with self.indent():
                for type_, event_list in by_pile_up_type.items():
                    if type_ == "gap":
                        self.li(f"skipped {len(event_list)} gap(s)")
                    else:
                        self.li(f"skipped {len(event_list)} {type_} extension(s)")

    def report_output_phase(self):
        tot_lost = total_len(
            self.events_by_reason,
            "maxInsertionError",
            "joinPolicy",
            "minExtensionLength",
        )
        self.li(f"lost {tot_lost} in `output` phase")

        with self.indent():
            self.report_max_insertion_error()
            self.report_join_policy()
            self.report_min_extension_length()

    def report_max_insertion_error(self):
        events = self.events_by_reason.pop("maxInsertionError", [])
        self.li(
            f"skipped {len(events)} insertion(s) because of {self.cli_opt('maxInsertionError')}"
        )
        by_insertion_type = group_by(events, lambda e: e.insertion.type)
        with self.indent():
            for type_, event_list in by_insertion_type.items():
                if type_ == "gap":
                    self.li(f"skipped {len(event_list)} gap(s)")
                else:
                    self.li(f"skipped {len(event_list)} {type_} extension(s)")

                if self.verbose:
                    with self.indent():
                        for event in event_list:
                            self.li(f"gap {self.gap(event.insertion)}")

    def report_join_policy(self):
        events = self.events_by_reason.pop("joinPolicy", [])
        self.li(
            f"skipped {len(events)} insertion(s) because of {self.cli_opt('joinPolicy')}"
        )
        by_insertion_type = group_by(events, lambda e: e.insertion.type)
        with self.indent():
            for type_, event_list in by_insertion_type.items():
                if type_ == "gap":
                    self.li(f"skipped {len(event_list)} gap(s)")
                else:
                    self.li(f"skipped {len(event_list)} {type_} extension(s)")

                if self.verbose:
                    with self.indent():
                        for event in event_list:
                            self.li(f"gap {self.gap(event.insertion)}")

    def report_min_extension_length(self):
        events = self.events_by_reason.pop("minExtensionLength", [])
        self.li(
            f"skipped {len(events)} extension(s) because of {self.cli_opt('minExtensionLength')}"
        )
        by_insertion_type = group_by(events, lambda e: e.insertion.type)
        with self.indent():
            for type_, event_list in by_insertion_type.items():
                self.li(f"skipped {len(event_list)} {type_} extension(s)")
                if self.verbose:
                    with self.indent():
                        for event in event_list:
                            self.li(f"{type_} extension of {self.gap(event.insertion)}")

    def report_unhandled_events(self):
        if len(self.events_by_reason) == 0:
            return

        self.li("unhandled skip-events:")
        with self.indent():
            for reason, event_list in self.events_by_reason.items():
                self.li(f"{reason} ({len(event_list)} times)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Search DENTIST logs for reasons why gaps are not closed."
    )
    parser.add_argument(
        "log_dir",
        metavar="<log-dir>",
        type=Path,
        help="DENTIST log directory (typically `./logs`)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="increase level of detail in the report",
    )

    return parser.parse_args()


def generate_report(log_dir, *, verbose=False, file=sys.stdout):
    report = Reporter(log_dir, verbose=verbose)
    report(file)


def main():
    args = parse_args()
    generate_report(args.log_dir, verbose=args.verbose)


if __name__ == "__main__":
    main()
    # TODO publish script on github
    # TODO reply to schraderl
