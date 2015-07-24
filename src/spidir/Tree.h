/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011
  
  Tree datastructure

=============================================================================*/

#ifndef SPIDIR_TREE_H
#define SPIDIR_TREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "ExtendArray.h"


using namespace std;


/*

    Parent Array Format (ptree)
    

           4
          / \
         3   \
        / \   \
       /   \   \
       0   1    2

     Then the parent tree array representation is
       ptree = [3, 3, 4, 4, -1]
     such that ptree[node's id] = node's parent's id

     In addition, the following must be true
     1. tree must be binary: n leaves, n-1 internal nodes
     2. leaves must be numbered 0 to n-1
     3. internal nodes are numbered n to 2n-2
     4. root must be numbered 2n-2
     5. the parent of root is -1
     6. the length of ptree is 2n-1

*/



namespace spidir {


// invariant: node->child->prev == last child of node
// invariant: [last child of node].next == NULL


// A node in the phylogenetic tree
class Node
{
public:
    Node(int nchildren=0) :
        name(-1),
        parent(NULL),
        //child(NULL),
        //next(NULL),
        //prev(NULL),
        children(NULL),
        nchildren(nchildren),
        dist(0.0)
    {
        if (nchildren != 0)
            allocChildren(nchildren);
    }
    
    ~Node()
    {
        if (children)
            delete [] children;
    }
    
    // Sets and allocates the number of children '_nchildren'
    void setChildren(int _nchildren)
    {
        children = resize(children, nchildren, _nchildren);
        nchildren = _nchildren;
    }
    
    // Allocates the number of children '_nchildren'
    void allocChildren(int _nchildren)
    {
        children = new Node* [_nchildren];
    }
    
    // Returns whether the node is a leaf
    bool isLeaf() const {
        return nchildren == 0;
    }
    
    // Adds a node 'node' to be a child
    void addChild(Node *node)
    {
        setChildren(nchildren + 1);
        children[nchildren - 1] = node;
        node->parent = this;

        /*
        // update pointers
        if (!child) {
            // adding first child
            child = node;
            node->prev = node;
        } else {
            Node *last = child->prev;
            last->next = node;
            node->prev = last;
            node->next = NULL;
            child->prev = node;
        }
        */
    }

    /*
    void RemoveChild(Node *node)
    {
        if (node->next) {
            // introduce your right sibling your left sibling
            nodem->next->prev = node->prev;
        } else if (elm != m_child) {
            // removing last child, update last pointer
            child->prev = node->prev;
        }
        
        if (node != child) {
            // introduce your left sibling your right sibling
            node->prev->next = node->next;
        } else {
            // if removing first child, find new first child
            // NOTE: node->next could be NULL and that is OK
            child = node->next;
        }
        
        // remove old links
        node->parent = NULL;
        node->next = NULL;
        node->prev = NULL;
    }
    */

    /*
    void ReplaceChild(Node *oldchild, Node *newchild)
    {
        // copy over links
        newchild->parent = this;
        newchild->next = oldchild->next;
        newchild->prev = oldchild->prev;
    
        // introduce newchild to siblings of old child
        if (newchild->next)
            newchild->next->prev = newchild;
        if (child == oldchild) {
            // replace first child
            child = newchild;

            if (oldchild->prev == oldchild)
                // replace single child
                newchild->prev = newchild;
        } else
            newchild->prev->next = newchild;
    
        oldchild->parent = NULL;
        oldchild->next = NULL;
        oldchild->prev = NULL;
    }
    */

    
    int name;           // node name id (matches index in tree.nodes)
    Node *parent;       // parent pointer
    //Node *child;        // first child
    //Node *next;         // next sibling
    //Node *prev;         // prev sibling

    Node **children;    // array of child pointers (size = nchildren)

    int nchildren;      // number of children
    float dist;         // branch length above node
    string longname;    // node name (used mainly for leaves only)
};


// A phylogenetic tree
class Tree
{
public:
    Tree(int nnodes=0) :
        nnodes(nnodes),
        root(NULL),
        nodes(nnodes, 100)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i] = new Node();
    }
    
    virtual ~Tree()
    {
        for (int i=0; i<nnodes; i++)
            delete nodes[i];
    }
    
    // Sets the branch lengths of the tree
    //  Arguments:
    //      dists: array of lengths (size = nnodes)
    void setDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i]->dist = dists[i];
    }
    
    
    // Gets the branch lengths of the tree
    //  Arguments:
    //      dists: output array (size = nnodes) for storing branch lengths
    void getDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            dists[i] = nodes[i]->dist;
    }
    
    // Sets the leaf names of the tree
    //  Arguments:
    //      names:      array of names (size > # leaves)
    //      leavesOnly: whether to only set names for leaves, or all nodes
    void setLeafNames(string *names, bool leavesOnly=true)
    {
        for (int i=0; i<nnodes; i++) {
            if (leavesOnly && !nodes[i]->isLeaf())
                nodes[i]->longname = "";
            else
                nodes[i]->longname = names[i];
        }
    }
    
    void reorderLeaves(string *names);
    
    // Gets leaf names of the nodes of a tree
    // Internal nodes are often named "" (empty string)
    //  Arguments:
    //      names:      output array for storing node names (size > nnodes)
    //      leavesOnly: whether to only get names for leaves, or all nodes
    void getNames(string *names, bool leavesOnly=true)
    {
        for (int i=0; i<nnodes; i++)
            if (!leavesOnly || nodes[i]->isLeaf())
                names[i] = nodes[i]->longname;
    }
        
    // Returns whether tree is rooted
    bool isRooted()
    {
        return (root != NULL && root->nchildren == 2);
    }
    
    // Returns the pointer to a node given is name id 'name'
    Node *getNode(int name)
    {
        return nodes[name];
    }
    
    // Adds a node 'node' to the tree
    // This will also set the node's name id
    Node *addNode(Node *node)
    {
        nodes.append(node);
        node->name = nodes.size() - 1;
        nnodes = nodes.size();
        return node;
    }
    
    // Compute a topology hash of the tree
    //  Arguments:
    //      key: output array (size = nnodes) containing a unique sequence of
    //           integers for the tree
    void hashkey(int *key);
    
    bool sameTopology(Tree *other);
    
    // set the topology to match another tree
    void setTopology(Tree *other);
    
    // Roots the tree on branch 'newroot'
    void reroot(Node *newroot, bool onBranch=true);
    
    // Roots the tree on branch connecting 'node1' and 'node2'
    void reroot(Node *node1, Node *node2);
    
    // Returns a new copy of the tree
    Tree *copy();
       
    // Returns whether the tree is self consistent
    bool assertTree();
        
public:    
    int nnodes;                 // number of nodes in tree
    Node *root;                 // root of the tree (NULL if no nodes)
    ExtendArray<Node*> nodes;   // array of nodes (size = nnodes)
};
                               

// A hash function for a topology key to an integer
struct HashTopology {
    static unsigned int hash(const ExtendArray<int> &key)
    {
        unsigned int h = 0, g;
        
        for (int i=0; i<key.size(); i++) {
            h = (h << 4) + key[i];
            if ((g = h & 0xF0000000))
                h ^= g >> 24;
            h &= ~g;
        }
        
        return h;
    }    
};


//=============================================================================
// Tree traversals

void getTreeSortedPostOrder(Tree *tree, ExtendArray<Node*> *nodes, 
                            int *ordering, Node *node=NULL);
void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);
void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);


//=============================================================================
// Input/output

void printFtree(int nnodes, int **ftree);
void printTree(Tree *tree, Node *node=NULL, int depth=0);


// C exports
extern "C" {


// Creates a 'forward tree' from a 'parent tree'
void makeFtree(int nnodes, int *ptree, int ***ftree);

// Deallocates a 'forward tree'
void freeFtree(int nnodes, int **ftree);


// Creates a tree object from a 'parent tree' array
void ptree2tree(int nnodes, int *ptree, Tree *tree);

// Creates a 'parent tree' from a tree object
void tree2ptree(Tree *tree, int *ptree);


Tree *makeTree(int nnodes, int *ptree);
void deleteTree(Tree *tree);
void setTreeDists(Tree *tree, float *dists);

}


} // namespace spidir


#endif // SPDIR_TREE_H

