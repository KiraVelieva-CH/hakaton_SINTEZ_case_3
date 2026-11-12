class GenomicDetectiveError(Exception):
    """Base application exception."""


class InvalidInputError(GenomicDetectiveError):
    """Input data are missing or malformed."""


class InvalidVCFError(InvalidInputError):
    """VCF structure or coordinates are invalid."""


class ModelExecutionError(GenomicDetectiveError):
    """An optional model or external tool failed."""


class ResourceUnavailableError(GenomicDetectiveError):
    """An optional annotation resource is unavailable."""
