#include "mesh.h"
///////////////////////////////////////////////////////////////////////////////

void dfs(const DualGraph *graph, Vertex *srcroot, deque<NodeType> &listQ) {
}

///////////////////////////////////////////////////////////////////////////////

void matchnodes(BinaryTree *tree, BinaryNode *child, BinaryNode *parent, 
                vector<NodePair> &matching) 
{
    if (!child->isMatched() && !parent->isMatched()) {
//      cout << " Match : " << child->getID() << " " << parent->getID() << endl; getchar();
        NodePair nodepair;
        nodepair.first = child;
        nodepair.second = parent;
        matching.push_back(nodepair);
        child->setMatchMark(1);
        parent->setMatchMark(1);
    }

    tree->removeNode(child);
    tree->removeNode(parent);

    return;
}

///////////////////////////////////////////////////////////////////////////////

void matchnode(BinaryTree *tree, BinaryNode* v, vector<NodePair> &matching) 
{
    BinaryNode *parv = v->getParent();

    if (parv == NULL) return;

    int degree = parv->getDegree();

    //Case 0: if par(v) is a node of degree 1, then T consists of
    //        two nodes joined by an edge. In this case, match v
    //        and par(v) and remove them from the tree. If par(v)
    //        is NIL then T consists of a single node and we leave
    //        it unmatched.

    if (degree == 1) {
        matchnodes(tree, v, parv, matching);
        return;
    }

    //Case 2: par(v) is a node of degree 2, match v and par(v)
    if (degree == 2) {
        matchnodes(tree, v, parv, matching);
        return;
    }

    //Case 3: par(v) is a node of degree 2 and v is the left
    //        child of par(v), In this case, match v and par(v)

    if ((degree == 3)) {
        BinaryNode *vsib = v->getSibling();
        assert(vsib);
        Vertex *d0 = v->getDualNode();
        Vertex *d1 = vsib->getDualNode();
        if (d0->getDegree() < d1->getDegree()) {
            matchnodes(tree, v, parv, matching);
            tree->removeNode(vsib);
        } else {
            matchnodes(tree, vsib, parv, matching);
            tree->removeNode(v);
        }
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////

bool has_same_dual(const BinaryNode *nd1, const BinaryNode *nd2) {
    return nd1->getID() < nd2->getID();
}
///////////////////////////////////////////////////////////////////////////////

bool already_matched(const BinaryNode *node) {
    return node->isMatched();
}
///////////////////////////////////////////////////////////////////////////////

BinaryNode* getChildofDegreeNParent(list<BinaryNode*> &levelnodes, int nd) 
{
    BinaryNode *currnode, *parent, *child;

    int ncount;

    BOOST_FOREACH(currnode, levelnodes) {
        parent = currnode->getParent();
        if (parent) {
            if (!parent->isMatched()) {
                ncount = 0;
                if (parent->getParent()) ncount = 1;
                for (int i = 0; i < parent->getNumChildren(); i++) {
                    child = parent->getChild(i);
                    if (!child->isMatched()) ncount++;
                }
                if (ncount == nd) return currnode;
            }
        }
    }

    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

BinaryNode *getNextNode(BinaryTree *tree, list<BinaryNode*> &levelnodes) 
{
    BinaryNode *currnode = NULL;

    if (levelnodes.empty()) return currnode;

    BOOST_FOREACH(currnode, levelnodes) {
        if (currnode->isMatched()) tree->removeNode(currnode);
    }


    list<BinaryNode*>::iterator it;
    it = remove_if(levelnodes.begin(), levelnodes.end(), already_matched);
    levelnodes.erase(it, levelnodes.end());

    BinaryNode *child0, *child1;

    // High Priority: parent having degree = 1;
    child0 = getChildofDegreeNParent(levelnodes, 1);

    if (!child0)
        child0 = getChildofDegreeNParent(levelnodes, 2);

    // Low Priority: parent having degree = 3;
    if (!child0)
        child0 = getChildofDegreeNParent(levelnodes, 3);

    return child0;
}

////////////////////////////////////////////////////////////////////////////////

void prunelevel(BinaryTree *tree, list<BinaryNode*> &levelnodes,
        vector<NodePair> &matching) 
{
    while (1) {
        BinaryNode *currnode = getNextNode(tree, levelnodes);
        if (currnode == NULL) break;
        matchnode(tree, currnode, matching);
    }
}

////////////////////////////////////////////////////////////////////////////////

void maximum_matching(BinaryTree *tree, vector<FacePair> &matching) 
{
    vector<NodePair> nodematching;

    matching.clear();
    int height = tree->getHeight();
    list<BinaryNode*> levelnodes;

    //Reset all the Matching marks to 0;
    for (int i = 0; i < height; i++) {
        levelnodes = tree->getLevelNodes(height - i - 1);
        BinaryNode *currnode;
        BOOST_FOREACH(currnode, levelnodes)
        currnode->setMatchMark(0);
    }

    for (int i = 0; i < height; i++) {
        levelnodes = tree->getLevelNodes(height - i - 1);
        prunelevel(tree, levelnodes, nodematching);
    }

    int nummatches = nodematching.size();
    matching.resize(nummatches);
    for (int i = 0; i < nummatches; i++) {
        matching[i].first = nodematching[i].first->getID();
        matching[i].second = nodematching[i].second->getID();
    }
}
///////////////////////////////////////////////////////////////////////////////

void maximum_matching(Mesh *inmesh, vector<FacePair> &matching) {

    cout << " Creating Dual Graph ... " << endl;
    DualGraph dgraph;
    dgraph.build(inmesh);

    inmesh->saveAs( "inmesh");

    cout << " Building Binary Tree of Dual Graph ... " << endl;
    BinaryTree *tree = new BinaryTree(&dgraph);
    tree->build();
    tree->saveAs("tree");

    cout << " Graph Matching ... " << endl;
    maximum_matching(tree, matching);

    cout << " maximum tree matching size ... " << matching.size() << endl;

    for (int i = 0; i < matching.size(); i++) {
        int id0 = matching[i].first;
        int id1 = matching[i].second;
        Vertex *du0 = inmesh->getFace(id0)->getDualNode();
        Vertex *du1 = inmesh->getFace(id1)->getDualNode();
        du0->setDualMate(du1);
        du1->setDualMate(du0);
    }

}
///////////////////////////////////////////////////////////////////////////////
