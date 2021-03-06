#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    LineGraph,
    LineGraphNode,
    LineGraphEdge,
    linegraph


struct LineGraphNode <: AbstractNode
    n1::Int
    n2::Int
end


struct LineGraphEdge <: UndirectedEdge
    u::Int
    v::Int
    node::Int
end

LineGraph = MapUDGraph{LineGraphNode, LineGraphEdge}


function linegraph(G::UDGraph)
    L = LineGraph()
    for (i, edge) in edgesiter(G)
        L.nodes[i] = LineGraphNode(edge.u, edge.v)
        L.adjacency[i] = Dict()
    end
    ecnt = 1
    for i in nodekeys(G)
        if degree(G, i) <= 1
            continue
        end
        # TODO: is there a stuff like itertools.combination?
        for e1 in neighboredgekeys(G, i)
            for e2 in neighboredgekeys(G, i)
                if e1 < e2
                    L.edges[ecnt] = LineGraphEdge(e1, e2, i)
                    L.adjacency[e1][e2] = ecnt
                    L.adjacency[e2][e1] = ecnt
                    ecnt += 1
                end
            end
        end
    end
    return L
end
