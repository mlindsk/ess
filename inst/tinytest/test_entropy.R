# Test that the two entropy implementations give same results
expect_equal(efs:::joint_entropy(cars), efs:::joint_entropy2(cars))
