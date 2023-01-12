import collections

class Queue:
    def __init__(self):
        self.elements = collections.deque()

    def empty(self):
        return len(self.elements) == 0

    def put(self, x):
        self.elements.append(x)

    def get(self):
        return self.elements.popleft()

class Graph:
    def __init__(self,bonds):
        self.__bonds = bonds

    @property
    def bonds(self):
        return self.__bonds

    def bfs(self, start):
        # печать того, что мы нашли
        frontier = Queue()
        frontier.put(start)
        visited = {}
        visited[start] = True
        points = []

        while not frontier.empty():
            current = frontier.get()
            points.append(current)
            for next in self.__bonds[current]:
                if next not in visited:
                    frontier.put(next)
                    visited[next] = True
        return points[1::]
