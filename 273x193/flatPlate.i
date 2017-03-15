Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-3
    max_iterations: 50
    kspace: 50
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

realms:

  - name: realm_1
    mesh: grid_struct_273x193_vol_ndtw.exo
    use_edges: yes
    check_for_missing_bcs: yes
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 100.0
     time_step_change_factor: 1.2

    equation_systems:
      name: theEqSys
      max_iterations: 2 

      solver_system_specification:
        velocity: solve_scalar
        turbulent_ke: solve_scalar
        specific_dissipation_rate: solve_scalar
        pressure: solve_cont

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-5

        - ShearStressTransport:
            name: mySST 
            max_iterations: 1
            convergence_tolerance: 1e-5

    initial_conditions:
      - constant: ic_1
        target_name: Unspecified-2-HEX
        value:
          pressure: 0
          velocity: [78.4197,0.0,0.0]
          turbulent_ke: 0.00138367
          specific_dissipation_rate: 9802.4625

    material_properties:
      target_name: Unspecified-2-HEX
      specifications:
        - name: density
          type: constant
          value: 1.177
        - name: viscosity
          type: constant
          value: 1.846e-5

    boundary_conditions:

    - wall_boundary_condition: bc_wall
      target_name: bottomwall
      wall_user_data:
        velocity: [0,0,0]
        turbulent_ke: 0.0
        use_wall_function: no

    - inflow_boundary_condition: bc_inflow
      target_name: inlet
      inflow_user_data:
        velocity: [78.4197,0.0,0.0]
        turbulent_ke: 0.00138367
        specific_dissipation_rate: 9802.4625

    - open_boundary_condition: bc_open
      target_name: outlet
      open_user_data:
        velocity: [0,0,0]
        pressure: 0.0
        turbulent_ke: 0.00138367
        specific_dissipation_rate: 9802.4625

    - symmetry_boundary_condition: bc_symBottom
      target_name: bottomsym
      symmetry_user_data:

    - open_boundary_condition: bc_open
      target_name: top
      open_user_data:
        velocity: [78.4197,0,0]
        pressure: 0.0
        turbulent_ke: 0.00138367
        specific_dissipation_rate: 9802.4625

    - periodic_boundary_condition: bc_front_back
      target_name: [front, back]
      periodic_user_data:
        search_tolerance: 0.0001

    solution_options:
      name: myOptions
      turbulence_model: sst

      options:
        - hybrid_factor:
            velocity: 1.0 
            turbulent_ke: 1.0
            specific_dissipation_rate: 1.0

        - alpha_upw:
            velocity: 1.0 

        - limiter:
            pressure: no
            velocity: yes
            turbulent_ke: yes
            specific_dissipation_rate: yes

        - projected_nodal_gradient:
            velocity: element
            pressure: element 
            turbulent_ke: element
            specific_dissipation_rate: element
    
        - input_variables_from_file:
            minimum_distance_to_wall: ndtw


    post_processing:
    
    - type: surface
      physics: surface_force_and_moment
      output_file_name: results/flatPlate.dat
      frequency: 25 
      parameters: [0,0]
      target_name: bottomwall

    output:
      output_data_base_name: results/flatPlate.e
      output_frequency: 100
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - pressure_force
       - tau_wall
       - turbulent_ke
       - specific_dissipation_rate
       - minimum_distance_to_wall
       - sst_f_one_blending
       - turbulent_viscosity

    restart:
      restart_data_base_name: restart/flatPlate.rst
      output_frequency: 2500
     
Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      time_step: 1.0e-10
      termination_time: 10.0 
      time_stepping_type: adaptive
      time_step_count: 0

      realms: 
        - realm_1
