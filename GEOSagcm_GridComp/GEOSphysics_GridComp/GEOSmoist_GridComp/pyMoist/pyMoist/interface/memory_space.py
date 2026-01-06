import enum


class MemorySpace(enum.Enum):
    CPU = 0
    GPU = 1


class InterfaceTransferType(enum.Enum):
    CPU_COPY = enum.auto()  # Copies because of layout mismatch
    CPU_MAP = enum.auto()  # No copy - memory map - same layout
    CPU_TO_GPU_TO_CPU = enum.auto()
