#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://https://github.com/watertap-org/reaktoro_enabled_watertap"
#################################################################################

import re
import sys

__author__ = "Alexander V. Dudchenko"


def get_ipopt_performance_data(solver_output_file):
    """process ipopt output file to extract performance data
    Args:
        solver_output_file (str): path to ipopt solver output file
    Returns:
        iters_match: regex match object with counts of various ipopt evals
        parsed_data: dict with parsed ipopt performance data"""
    with open(solver_output_file) as f:
        lines = f.read()
    iters_match = re.search(
        r"""
        Number\ of\ objective\ function\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ objective\ gradient\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ equality\ constraint\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ inequality\ constraint\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ equality\ constraint\ Jacobian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ inequality\ constraint\ Jacobian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        Number\ of\ Lagrangian\ Hessian\ evaluations\s*=\s*([-+eE0-9.]+)\s*
        """,
        lines,
        re.DOTALL | re.VERBOSE,
    )
    return iters_match, _parse_ipopt_output(lines)


def _parse_ipopt_output(output):
    _ALPHA_PR_CHARS = set("fFhHkKnNRwstTr")

    parsed_data = {}
    import io

    # Convert output to a string so we can parse it
    if isinstance(output, io.StringIO):
        output = output.getvalue()

    # Extract number of iterations
    iter_match = re.search(r"Number of Iterations.*:\s+(\d+)", output)
    if iter_match:
        parsed_data["iters"] = int(iter_match.group(1))
    # Gather all the iteration data
    iter_table = re.findall(r"^(?:\s*\d+.*?)$", output, re.MULTILINE)
    if iter_table:
        columns = [
            "iter",
            "objective",
            "inf_pr",
            "inf_du",
            "lg_mu",
            "d_norm",
            "lg_rg",
            "alpha_du",
            "alpha_pr",
            "ls",
        ]
        iterations = []

        for line in iter_table:
            tokens = line.strip().split()
            if len(tokens) != len(columns):
                continue
            iter_data = dict(zip(columns, tokens))

            # Extract restoration flag from 'iter'
            iter_data["restoration"] = iter_data["iter"].endswith("r")
            if iter_data["restoration"]:
                iter_data["iter"] = iter_data["iter"][:-1]

            # Separate alpha_pr into numeric part and optional tag
            iter_data["step_acceptance"] = iter_data["alpha_pr"][-1]
            if iter_data["step_acceptance"] in _ALPHA_PR_CHARS:
                iter_data["alpha_pr"] = iter_data["alpha_pr"][:-1]
            else:
                iter_data["step_acceptance"] = None

            # Attempt to cast all values to float where possible
            for key in columns:
                if iter_data[key] == "-":
                    iter_data[key] = None
                else:
                    try:
                        iter_data[key] = float(iter_data[key])
                    except (ValueError, TypeError):
                        print(
                            "Error converting Ipopt log entry to "
                            f"float:\n\t{sys.exc_info()[1]}\n\t{line}"
                        )

            # assert len(iterations) == iter_data.pop("iter"), (
            #     f"Parsed row in the iterations table\n\t{line}\ndoes not "
            #     f"match the next expected iteration number ({len(iterations)})"
            # )
            iterations.append(iter_data)

        parsed_data["iteration_log"] = iterations

    # Extract scaled and unscaled table
    scaled_unscaled_match = re.search(
        r"""
        Objective\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Dual\ infeasibility\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Constraint\ violation\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        (?:Variable\ bound\ violation:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*)?
        Complementarity\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)\s*
        Overall\ NLP\ error\.*:\s*([-+eE0-9.]+)\s+([-+eE0-9.]+)
        """,
        output,
        re.DOTALL | re.VERBOSE,
    )

    if scaled_unscaled_match:
        groups = scaled_unscaled_match.groups()
        all_fields = [
            "incumbent_objective",
            "dual_infeasibility",
            "constraint_violation",
            "variable_bound_violation",  # optional
            "complementarity_error",
            "overall_nlp_error",
        ]

        # Filter out None values and create final fields and values.
        # Nones occur in old-style IPOPT output (<= 3.13)
        zipped = [
            (field, scaled, unscaled)
            for field, scaled, unscaled in zip(all_fields, groups[0::2], groups[1::2])
            if scaled is not None and unscaled is not None
        ]

        scaled = {k: float(s) for k, s, _ in zipped}
        unscaled = {k: float(u) for k, _, u in zipped}

        parsed_data.update(unscaled)
        parsed_data["final_scaled_results"] = scaled

    # Newer versions of IPOPT no longer separate timing into
    # two different values. This is so we have compatibility with
    # both new and old versions
    parsed_data["cpu_seconds"] = {
        k.strip(): float(v)
        for k, v in re.findall(
            r"Total(?: CPU)? sec(?:ond)?s in ([^=]+)=\s*([0-9.]+)", output
        )
    }

    return parsed_data
