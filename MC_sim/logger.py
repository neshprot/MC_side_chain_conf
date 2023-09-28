"""
create log file
"""
from abc import abstractmethod, ABC


class Logger(ABC):
    """
    abstract class for log
    """
    def __init__(self):
        pass

    def __call__(self, *msg):
        self.write(*msg)

    @abstractmethod
    def write(self, *msg):
        """
        writing smth into log file
        :param msg:
        :return:
        """


class ConsoleLogger(Logger):
    """
    print into console
    """
    def write(self, *msg):
        print(*msg)


class FileLogger(Logger):
    """
    class for log file
    """
    def __init__(self, fname):
        super().__init__()
        self.fname = fname

        with open(fname, "w", encoding="utf-8") as file:
            file.write("\n")

    def write(self, *msg):
        with open(self.fname, "a", encoding="utf-8") as file:
            print(*msg, file=file, end="", sep="")
