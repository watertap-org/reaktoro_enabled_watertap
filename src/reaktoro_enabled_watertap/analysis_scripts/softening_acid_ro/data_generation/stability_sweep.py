from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
import reaktoro_enabled_watertap.flowsheets.softening_acid_ro.softening_acid_ro as sar


def solve_with_ma27(m, tee=True, **kwargs):
    result = sar.solve_model(m, tee=tee, linear_solver="ma27")
    return result


def main():
    cwd = get_working_dir()
    loopTool(
        cwd + "/stability_sweep.yaml",
        build_function=sar.build_model,
        initialize_function=sar.initialize,
        optimize_function=solve_with_ma27,
        save_name="stability_sweep",
        probe_function=sar.test_func,
        saving_dir=cwd,
        number_of_subprocesses=1,
        num_loop_workers=3,
    )


if __name__ == "__main__":
    main()
