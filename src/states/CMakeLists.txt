add_fsm_state_simple(FSMSequenceOne_Initial)

add_fsm_state_simple(Operation)
target_link_libraries(Operation PUBLIC mc_iam::devices)

add_fsm_data_directory(data)
