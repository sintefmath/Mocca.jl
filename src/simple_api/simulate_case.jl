function simulate_case(case)

    return Jutul.simulate!(case.sim, case.timesteps;
    config = case.cfg,
    forces = case.sim_forces,
);

end