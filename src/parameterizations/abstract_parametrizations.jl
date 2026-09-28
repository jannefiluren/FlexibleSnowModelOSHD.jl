abstract type AbstractParameterization{Tf <: Real} end

abstract type AbstractAlbedo{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractCompaction{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractConductivity{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractFreshSnowDensity{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractHydrology{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractLayering{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractSnowFraction{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractStabilityCorrection{Tf} <: AbstractParameterization{Tf} end

abstract type AbstractLandCover{Tf} <: AbstractParameterization{Tf} end
abstract type AbstractSubstrate{Tf} <: AbstractParameterization{Tf} end
