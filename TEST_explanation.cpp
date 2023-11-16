// (0.3) Check that the simulation produces the mathematically correct answer when u = 0 and v = 0.

中文解释：
概述：检查再给定u、v固定值0的情况下是否能输出正确结果
步骤：
在给定初始条件下，如果u和v都被初始化为0，并且在simulateStep()函数中的迭代过程中不发生任何变化，那么在每次迭代后u和v的值将始终保持为0。
这是因为在计算下一个时间步长时，dU和dV的计算都依赖于当前的u和v值。由于u和v始终为0，那么dU和dV也将始终为0，因此下一个时间步长的u和v值仍然为0。
另外，根据countElementsAboveThreshold()函数的定义，该函数用于计算v中大于阈值的元素的比例。由于在给定的情况下v始终为0，因此不存在大于阈值的元素，所以返回的比例将始终为0。

TEST(simulation, check_simulation_with_zero_initial_condition) {
    // Set initial conditions to u = 0 and v = 0
    u = std::vector<std::vector<double>>(width, std::vector<double>(height, 0.0));
    v = std::vector<std::vector<double>>(width, std::vector<double>(height, 0.0));
    F=0.0;
    k=0.0;
    for (int iteration = 0; iteration < numIterations; ++iteration) {
        simulateStep();
    }
    // check all zero
    for (int i = 0; i < u.size(); i++)
    {
        for (int j = 0; j < u[i].size(); j++)
        {
            ASSERT_EQ(u[i][j], 0);
        }
    }

    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            ASSERT_EQ(v[i][j], 0);
        }
    }
    // use countElementsAboveThreshold 
    double proportionAboveThreshold = countElementsAboveThreshold(threshold);

    ASSERT_EQ(proportionAboveThreshold, 0.0);
}

// # additional tests:

// (0.4) Check initial values of u and v

中文解释：基于初始函数，检验是否：u 在任何点都是1，v 只在固定的区域内为随机值，且这个随机值在 [0,1) 间

TEST(simulation, check_initial_values) {
    // Run initialization
    init();

    // Check that u has a value of 1.0 and v is within (0, 1) range for the specified region
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            ASSERT_EQ(u[x][y],1);
            if (x >= 100 && x <= 200 && y >= 100 && y <= 150) {
                ASSERT_GE(v[x][y], 0.0);
                ASSERT_LT(v[x][y], 1.0);
            } else {
                ASSERT_EQ(v[x][y], 0.0);
            }
        }
    }
}

// (0.5) Check if simulation step is working correctly and u and v are updated as expected

中文解释：
概述：simulation函数是否正常运行
步骤：首先，运行 init() 函数进行初始化操作，然后保存了 u 和 v 的初始状态。
接下来，调用 simulateStep() 函数执行一次模拟的步骤。
之后，通过比较 u 和 v 的当前状态与初始状态的值，检查它们是否有更新。如果在 u 和 v 的任何一个元素上发现了变化（不等于初始状态的对应元素），则将 isUpdated 标志设置为 true。
最后，使用断言 ASSERT_TRUE(isUpdated) 来验证 isUpdated 是否为 true。如果 isUpdated 为 true，表示 u 和 v 已经被成功更新，测试通过。如果 isUpdated 为 false，表示 u 和 v 没有被更新，测试失败。

TEST(simulation, check_simulation_step) {
    // Run initialization
    init();

    // Save the initial state of u and v
    std::vector<std::vector<double>> initialU = u;
    std::vector<std::vector<double>> initialV = v;

    // Perform one simulation step
    simulateStep();

    // Check that u and v have been updated
    bool isUpdated = false;
    for (int x = 1; x < width - 1; ++x) {
        for (int y = 1; y < height - 1; ++y) {
            if (u[x][y] != initialU[x][y] || v[x][y] != initialV[x][y]) {
                isUpdated = true;
                break;
            }
        }
        if (isUpdated) {
            break;
        }
    }
    ASSERT_TRUE(isUpdated);
}


// (0.6) Test countElementsAboveThreshold function

中文解释：
概述：用于检查 countElementsAboveThreshold() 函数的功能是否正确。

首先，通过将 v 初始化为一个大小为 width * height 的二维向量，并将所有元素的值设置为0.0。
然后，为了测试 countElementsAboveThreshold 函数，将特定的值（0.2）赋给 v 的一部分元素。具体而言，将 v[x][y] 的值设置为0.2，其中 x 和 y 的取值范围是 0 到 100。
接下来，计算在 v 中大于阈值的元素个数 countAboveThreshold，并计算了总元素个数 totalElements（即 width * height）。
然后，使用 countElementsAboveThreshold 函数获取实际的比例 actualRatio。
接着，计算了预期的比例 expectedRatio，即 countAboveThreshold 除以 totalElements。
最后，使用断言 EXPECT_NEAR(actualRatio, expectedRatio, epsilon) 来验证实际比例 actualRatio 是否与预期比例 expectedRatio 在 epsilon 的范围内接近。epsilon 是一个小的误差范围，用于处理浮点数比较时的舍入误差。


TEST(simulation, check_count_elements_above_threshold) {
    // Set para
    v = std::vector<std::vector<double>>(width, std::vector<double>(height, 0.0));

    // Set a specific value for v[x][y] to test countElementsAboveThreshold function
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (x >= 0 && x <= 100 && y >= 0 && y <= 100){
                v[x][y] = 0.2;
            }else{
                v[x][y] = 0;
            }
        }
    }

    // Caulate the expected ratio
    int countAboveThreshold = 0;
    int totalElements = width * height;
    for (const auto& row : v) {
        for (double value : row) {
            if (value > threshold) {
                countAboveThreshold++;
            }
        }
    }
    double expectedRatio = static_cast<double>(countAboveThreshold) / totalElements;

    // Get the actual Ratio
    double actualRatio = countElementsAboveThreshold(threshold);

    // Compare the actual ratio with the expected ratio in epsilon
    double epsilon = 1e-6;
    EXPECT_NEAR(actualRatio, expectedRatio, epsilon);
}



// (0.7) Check if simulation correctly outputs VTK files

中文解释：概述：验证模拟过程是否正确地输出了 VTK 文件

具体步骤：
在每个输出间隔（outputInterval）的迭代步骤中，调用 writeVTKFile(iteration) 函数来将当前模拟状态写入 VTK 文件。
然后，检查生成的 VTK 文件是否存在。它通过构造文件名 "output_iteration.vtk"，并尝试打开该文件来进行检查。
如果文件存在（即 vtkFile.good() 返回 true），则断言 ASSERT_TRUE(vtkFile.good()) 成立，表示生成的 VTK 文件存在。如果文件不存在或打开失败，则断言失败，测试失败。

TEST(simulation, check_vtk_output) {
    // Run initialization
    init();

    // Perform simulations
    for (int iteration = 0; iteration < numIterations; ++iteration) {
        simulateStep();

        // Write VTK file for each output interval
        if (iteration % outputInterval == 0) {
            writeVTKFile(iteration);

            // Check that the VTK file exists
            std::stringstream ss;
            ss << "output_" << iteration << ".vtk";
            std::ifstream vtkFile(ss.str());
            ASSERT_TRUE(vtkFile.good());
            vtkFile.close();
        }
    }
}