//
// Created by Truong Giang Do on 02/12/2023.
//

#ifndef SSSP_NEW_CUSTOMPRIORITYQUEUE_H
#define SSSP_NEW_CUSTOMPRIORITYQUEUE_H

#include "Graph.h"

template<typename Node>
class custom_priority_queue : public std::priority_queue<Node, std::vector<Node>> {
public:
    custom_priority_queue(size_t reserve_size) {
        this->c.reserve(reserve_size);
    }

    custom_priority_queue() {
    }
};

#endif //SSSP_NEW_CUSTOMPRIORITYQUEUE_H
