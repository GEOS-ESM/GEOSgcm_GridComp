from ndsl import OptimizationConfig, StencilFactory


def get_optimization_config(stencil_factory: StencilFactory) -> OptimizationConfig:
    return OptimizationConfig(
        stree=OptimizationConfig.Tree(
            enabled=True,
            merger=OptimizationConfig.Tree.Merger(enabled=True, overcompute=True),
            refine_transients=False if stencil_factory.backend.is_gpu_backend() else True,
        ),
    )
