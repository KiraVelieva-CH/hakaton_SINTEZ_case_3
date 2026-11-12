from abc import ABC, abstractmethod

class BaseAnalyzer(ABC):
    name = "base"
    @abstractmethod
    def analyze(self, context: dict) -> dict:
        raise NotImplementedError
