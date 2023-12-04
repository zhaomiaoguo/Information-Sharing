import math
from typing import List

class Node():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def copy(self):
        return Node(self.x, self.y)

class Edge():
    def __init__(self, node1, node2, length):
        self.node1 = node1
        self.node2 = node2
        self.length = length

class Network():

    def __init__(self, nodes: List[Node], edges: List[Edge]):
        self.nodes = nodes
        self.edges = edges

    def scale_network(self, scale):

        def scale_node(pivot: Node, node: Node) -> Node:
            
            base_len = math.dist((pivot.x,pivot.y), (node.x,node.y))
            new_len = base_len * scale

            theta = math.atan2(node.y-pivot.y, node.x-pivot.x)

            print(f'angle: {math.degrees(theta)}')

            x_dist = new_len * math.cos(theta)
            y_dist = new_len * math.sin(theta)

            return Node(pivot.x+x_dist,pivot.y+y_dist)


        pivot = self.nodes[0]

        scaled_nodes = [pivot.copy()]
        scaled_edges = []

        for node in self.nodes:
            if node == pivot: continue

            scaled_nodes.append(scale_node(pivot, node))

        for edge in self.edges:
            continue

        self.nodes = scaled_nodes

def main():
    nodes = [
        Node(-27,50),
        Node(13,50),
        Node(36.81,68.18),
        Node(98.88,26.80),
        Node(125.80, 51.54),
        Node(146.94, 51.72)
    ]

    # nodes = [Node(-15,32), Node(5,52), Node(5,32), Node(-15,52)]

    net = Network(nodes, edges=[])
    net.scale_network(scale=107.43)
    
    for node in net.nodes:
        print(f'x: {"%.2f"%node.x}, y: {"%.2f"%node.y}')

if __name__ == '__main__':
    main()