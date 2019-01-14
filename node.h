#ifndef NODE_H
#define NODE_H

#include <QPointF>

struct node{
    uint32_t id;
    QPointF coord;
    node() = default;
    node(uint32_t nid,QPointF ncoord):id(nid),coord(ncoord){}
};


#endif // NODE_H
