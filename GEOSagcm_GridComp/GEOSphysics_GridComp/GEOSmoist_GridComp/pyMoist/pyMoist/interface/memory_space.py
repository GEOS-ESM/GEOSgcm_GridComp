import enum


class MemorySpace(enum.Enum):
    CPU = 0
    GPU = 1


class InterfaceTransferType(enum.Enum):
    ALL_CPU = enum.auto()
    CPU_TO_GPU_TO_CPU = enum.auto()
