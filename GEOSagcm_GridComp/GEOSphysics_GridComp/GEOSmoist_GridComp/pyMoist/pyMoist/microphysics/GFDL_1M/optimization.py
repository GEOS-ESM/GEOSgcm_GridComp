from ndsl import OptimizationConfig, StencilFactory


def get_optimization_config(stencil_factory: StencilFactory):
    return OptimizationConfig(stree=OptimizationConfig.Tree(enabled=True, merger=OptimizationConfig.Tree.Merger(enabled=True, overcompute=True)))
