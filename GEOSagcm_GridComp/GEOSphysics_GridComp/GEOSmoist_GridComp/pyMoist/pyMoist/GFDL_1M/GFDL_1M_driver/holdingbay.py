            print("TESTING: precip_fall_ice")
        tot_indicies = 0
        failed_indicies = 0
        testing_variable = stencil.precip_fall
        reference_variable = precip_fall_ice
        for i in range(testing_variable.view[:].shape[0]):
            for j in range(testing_variable.view[:].shape[1]):
                tot_indicies = tot_indicies + 1
                if testing_variable.view[i, j, 0] != reference_variable.view[i, j]:
                    failed_indicies = failed_indicies + 1
                    print(
                        "precip_fall_ice DIFF: ",
                        i,
                        j,
                        "computed: ",
                        testing_variable.view[i, j, 0],
                        max(stencil.qi1.view[i, j, :]),
                        "reference: ",
                        reference_variable.view[i, j],
                        max(qi_terminal_fall.view[i, j, :]),
                        "difference: ",
                        testing_variable.view[i, j, 0] - reference_variable.view[i, j],
                    )
        print("precip_fall_ice failures: ", failed_indicies, "/", tot_indicies)

        print("TESING: qi")
        tot_indicies = 0
        failed_indicies = 0
        testing_variable = stencil.qi1
        reference_variable = qi_terminal_fall
        for i in range(testing_variable.view[:].shape[0]):
            for j in range(testing_variable.view[:].shape[1]):
                for k in range(testing_variable.view[:].shape[2]):
                    tot_indicies = tot_indicies + 1
                    if (
                        testing_variable.view[i, j, k]
                        != reference_variable.view[i, j, k]
                    ):
                        failed_indicies = failed_indicies + 1
                        print(
                            "DIFF: ",
                            i,
                            j,
                            k,
                            "computed: ",
                            testing_variable.view[i, j, k],
                            "reference: ",
                            reference_variable.view[i, j, k],
                            "difference: ",
                            testing_variable.view[i, j, k]
                            - reference_variable.view[i, j, k],
                        )
        print("failures: ", failed_indicies, "/", tot_indicies)