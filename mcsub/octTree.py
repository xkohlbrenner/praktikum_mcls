from node import Node

class OctTree:

    def __init__(self):
        self.Tree= []

    #creates octan tree and applies the environment to all nodes at the lowest area
    def createTree(self, x1, x2, y1, y2, z1, z2, depth, envList):
        tempbranch = Node(x1, x2, y1, y2, z1, z2)
        self.Tree.append(tempbranch)
        if depth == 0:
            tempbranch.set_env(envList[int(abs(x1)*64*64+abs(y1)*64+abs(z1)*64)])
        else:
            depth -= 1
            c1 = self.createTree(x1, x2-(x2-x1)/2, y1, y2-(y2-y1)/2, z1, z2-(z2-z1)/2, depth, envList)
            c2 = self.createTree(x1, x2-(x2-x1)/2, y1, y2-(y2-y1)/2, z2-(z2-z1)/2, z2, depth, envList)
            c3 = self.createTree(x1, x2-(x2-x1)/2, y2-(y2-y1)/2, y2, z1, z2-(z2-z1)/2, depth, envList)
            c4 = self.createTree(x1, x2-(x2-x1)/2, y2-(y2-y1)/2, y2, z2-(z2-z1)/2, z2, depth, envList)
            c5 = self.createTree(x2-(x2-x1)/2, x2, y1, y2-(y2-y1)/2, z1, z2-(z2-z1)/2, depth, envList)
            c6 = self.createTree(x2-(x2-x1)/2, x2, y1, y2-(y2-y1)/2, z2-(z2-z1)/2, z2, depth, envList)
            c7 = self.createTree(x2-(x2-x1)/2, x2, y2-(y2-y1)/2, y2, z1, z2-(z2-z1)/2, depth, envList)
            c8 = self.createTree(x2-(x2-x1)/2, x2, y2-(y2-y1)/2, y2, z2-(z2-z1)/2, z2, depth, envList)
            tempbranch.set_child(c1)
            tempbranch.set_child(c2)
            tempbranch.set_child(c3)
            tempbranch.set_child(c4)
            tempbranch.set_child(c5)
            tempbranch.set_child(c6)
            tempbranch.set_child(c7)
            tempbranch.set_child(c8)
            #look out if all children of the node have the same environment
            #if this is the case set the env of this node to same as of the children
            if all(child.env == tempbranch.children[0].env for child in tempbranch.children):
                tempbranch.env = tempbranch.children[0].env
            else:
                tempbranch.env = 0
        return tempbranch
    
    #a nodes gets an environment, if all children are part of the same environment
    #this process can first be done, if the tree is build and all nodes at the lowest area have a environment
    def set_environments(self):
        self.Tree[0].look_env()
    
    #returns the environment based on the porvided point variabless
    def get_env(self, x, y, z):
        #Tree[0] is the first node and therefore the root of the tree
        node = self.Tree[0]
        env = node.get_env(x, y, z)

        return env
        
    
    def get_node(self, node):
         return self.Tree[node]