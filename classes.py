setNumber = 0
edgeNumber = 0
class Set:
    def __init__(self, points, intersection):
        self.points = tuple(points)
        self.intersection = intersection
        self.algo2numhascrossed = 0
        self.algo3weight = 1
        global setNumber
        self.id = setNumber
        setNumber += 1

    @classmethod
    def frompointlist(cls,points):
        return cls(points,0)

    def increase(self):
        self.intersection += 1

    def decrease(self):
        self.intersection -= 1
    
    def algo2inc(self):
        self.algo2numhascrossed += 1

    def algo3inc(self):
        self.algo3weight *= 2

class Edge:
    def __init__(self, points, intersection):
        if len(points) != 2:
            raise ValueError('Edge can only have 2 pts, given:', len(points))
        self.points = tuple(points)
        self.intersection = intersection
        self.algo3weight = 0
        global edgeNumber
        self.id = edgeNumber
        edgeNumber += 1

    @classmethod
    def frompointlist(cls,points):
        if len(points) != 2:
            raise ValueError('Edge can only have 2 pts, given:', len(points))
        return cls(points,0)

    def increase(self):
        self.intersection += 1

    def decrease(self):
        self.intersection -= 1
