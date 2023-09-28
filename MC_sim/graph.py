"""
tools for creating graph
"""
import collections


class Queue:
    """
    queue for dots
    """
    def __init__(self):
        self.elements = collections.deque()

    def empty(self):
        """
        check if queue is empty
        :return: true or false
        """
        return len(self.elements) == 0

    def put(self, value):
        """
        put element in queue
        :param value:
        :return:
        """
        self.elements.append(value)

    def get(self):
        """
        get the element from queue
        :return:
        """
        return self.elements.popleft()


class Graph:
    """
    graph and tools for working with it
    """
    def __init__(self, bonds):
        self.__bonds = bonds

    @property
    def bonds(self):
        """
        get bonds
        :return:
        """
        return self.__bonds

    def bfs(self, start):
        """
        graph bfs
        :param start:
        :return:
        """
        frontier = Queue()
        frontier.put(start)
        visited = {}
        visited[start] = True
        points = []

        while not frontier.empty():
            current = frontier.get()
            points.append(current)
            for i in self.__bonds[current]:
                if i not in visited:
                    frontier.put(i)
                    visited[i] = True
        return points[1::]
