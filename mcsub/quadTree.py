from node import Node

class QuadTree:

    def __init__(self):
        self.Tree= []

    #creates octan tree and applies the environment to all nodes at the lowest area
    def create_Tree(self, r1, r2, h1, h2, pixel, depth, envList):
        tempbranch = Node(r1, r2, h1, h2)
        self.Tree.append(tempbranch)
        if depth == 0:
            #print(str(x1) + " " + str(h1) + " " + str(z1))
            #print(int(abs(x1)*pixel*pixel+abs(h1)*pixel+abs(z1)))
            half = int(pixel/2)
            tempbranch.set_env(envList[int((h1+half)*pixel+(r1+half))])
        else:
            depth -= 1
            c1 = self.create_Tree(r1, r2-(r2-r1)/2, h1, h2-(h2-h1)/2, pixel, depth, envList)
            c2 = self.create_Tree(r1, r2-(r2-r1)/2, h2-(h2-h1)/2, h2, pixel, depth, envList)
            c3 = self.create_Tree(r2-(r2-r1)/2, r2, h1, h2-(h2-h1)/2, pixel, depth, envList)
            c4 = self.create_Tree(r2-(r2-r1)/2, r2, h2-(h2-h1)/2, h2, pixel, depth, envList)
            tempbranch.set_child(c1)
            tempbranch.set_child(c2)
            tempbranch.set_child(c3)
            tempbranch.set_child(c4)
            #look out if all children of the node have the same environment
            #if this is the case set the env of this node to same as of the children
            if all(child.env == tempbranch.children[0].env for child in tempbranch.children):
                tempbranch.env = tempbranch.children[0].env
            else:
                tempbranch.env = 0
        return tempbranch
    
       
    #returns the environment based on the porvided point variabless
    def get_env(self, x, y):
        #Tree[0] is the first node and therefore the root of the tree
        node = self.Tree[0]
        env = node.get_env(x, y)

        return env
        
    
    def get_node(self, node):
         return self.Tree[node]