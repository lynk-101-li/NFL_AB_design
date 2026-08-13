"""Validated planners for real external antibody-design backends.

The adapters compile inputs and exact command specifications; they never run a
model implicitly. This keeps a hard boundary between simulated repository
outputs and real RFantibody, IgGM, or Germinal results.
"""

from .germinal import GerminalAdapterError, build_germinal_jobs
from .iggm import IgGMAdapterError, build_iggm_jobs
from .rfantibody import RFantibodyAdapterError, build_rfantibody_plan

__all__ = [
    "GerminalAdapterError",
    "IgGMAdapterError",
    "RFantibodyAdapterError",
    "build_germinal_jobs",
    "build_iggm_jobs",
    "build_rfantibody_plan",
]
