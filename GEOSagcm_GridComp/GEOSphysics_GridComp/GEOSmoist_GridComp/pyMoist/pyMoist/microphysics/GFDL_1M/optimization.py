from ndsl import OptimizationConfig


def get_optimization_config() -> OptimizationConfig:
    return OptimizationConfig(
        stree=OptimizationConfig.Tree(
            enabled=True,
            merger=OptimizationConfig.Tree.Merger(enabled=True, overcompute=True),
        ),
    )
