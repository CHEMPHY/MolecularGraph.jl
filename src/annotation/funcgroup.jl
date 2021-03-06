#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    funcgrouptable,
    FUNC_GROUP_TABLE,
    FGTermNode,
    FGRelationEdge,
    FunctionalGroup,
    functionalgroup!,
    fgrouprecord,
    fgroupcond,
    fgroupquery,
    largestcomponents


function funcgrouptable()
    table = []
    files = [
        "funcgroup.yaml",
        "ring.yaml",
        "biomolecule.yaml"
    ]
    dir = joinpath(dirname(@__FILE__), "..", "..", "assets", "funcgroup")
    for f in files
        src = joinpath(dir, f)
        data = YAML.load(open(src))
        println("loading: $(f)")
        append!(table, data)
    end
    return table
end

FUNC_GROUP_TABLE = funcgrouptable()


struct FGTermNode <: AbstractNode
    term::Symbol
end


struct FGRelationEdge <: DirectedEdge
    source::Int
    target::Int
    relation::Symbol
end


struct FunctionalGroup <: Annotation
    nodeset::Dict{Symbol,Set{Set{Int}}}
    graph::MapDGraph{FGTermNode,FGRelationEdge}

    function FunctionalGroup()
        new(Dict(), MapDGraph{FGTermNode,FGRelationEdge}())
    end
end


function functionalgroup!(mol::VectorMol; recalculate=false)
    if haskey(mol.annotation, :FunctionalGroup) && !recalculate
        return
    end
    aromatic!(mol, recalculate=recalculate)
    mol.annotation[:FunctionalGroup] = fg = FunctionalGroup()
    fggraphidx = Dict{Symbol,Int}()
    ncnt = 0
    ecnt = 0
    for rcd in FUNC_GROUP_TABLE
        fgset = fgrouprecord(mol, rcd)
        fgkey = Symbol(rcd["key"])
        fg.nodeset[fgkey] = fgset
        # Update ontology graph
        if !isempty(fgset)
            ncnt += 1
            fggraphidx[fgkey] = ncnt
            updatenode!(fg.graph, FGTermNode(fgkey), ncnt)
            if "have" in keys(rcd)
                for k in rcd["have"]
                    ecnt += 1
                    e = FGRelationEdge(fggraphidx[Symbol(k)], ncnt, :partof)
                    updateedge!(fg.graph, e, ecnt)
                end
            end
            if "isa" in keys(rcd)
                for k in rcd["isa"]
                    ecnt += 1
                    e = FGRelationEdge(ncnt, fggraphidx[Symbol(k)], :isa)
                    updateedge!(fg.graph, e, ecnt)
                end
            end
        end
    end
end


function fgrouprecord(mol::VectorMol, rcd)
    fgsetmap = mol.annotation[:FunctionalGroup].nodeset
    newset = Set{Set{Int}}()
    # Membership filter
    if "have" in keys(rcd)
        for k in rcd["have"]
            if isempty(fgsetmap[Symbol(k)])
                return newset
            end
        end
    end
    # Substructure match
    if "isa" in keys(rcd)
        q = parse(ConnectedSMARTS, rcd["query"])
        for k in rcd["isa"]
            refset = fgsetmap[Symbol(k)]
            eachset = Set{Set{Int}}()
            for s in refset
                subst = atomsubstr(mol, s)
                for (emap, nmap) in querymatchiter(subst, q, mode=:graph)
                    if !isempty(emap)
                        esub = edgesubgraph(mol.graph, keys(emap))
                        push!(eachset, nodekeys(esub))
                    elseif !isempty(nmap)
                        push!(eachset, keys(nmap))
                    end
                end
            end
            if isempty(newset)
                union!(newset, eachset)
            else
                intersect!(newset, eachset)
            end
        end
    else
        union!(newset, fgroupcond(mol, rcd))
    end
    return newset
end


function fgroupcond(mol::VectorMol, rcd)
    newset = Set{Set{Int}}()
    if "any" in keys(rcd)
        for query in rcd["any"]
            union!(newset, fgroupquery(mol, query))
        end
    else
        union!(newset, fgroupquery(mol, rcd["query"]))
    end
    return newset
end


function fgroupquery(mol::VectorMol, query)
    q = parse(ConnectedSMARTS, query)
    newset = Set{Set{Int}}()
    for (emap, nmap) in querymatchiter(mol, q)
        if !isempty(emap)
            esub = edgesubgraph(mol.graph, keys(emap))
            push!(newset, nodekeys(esub))
        elseif !isempty(nmap)
            push!(newset, keys(nmap))
        end
    end
    return newset
end



function largestcomponents(fg::FunctionalGroup)
    components = Dict{Symbol,Set{Set{Int}}}()
    ontnodes = topologicalsort(fg.graph)
    for n in ontnodes
        rmset = Set{Set{Int}}()
        nterm = getnode(fg.graph, n).term
        for (s, e) in successors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :partof
                continue
            end
            ansterm = getnode(fg.graph, s).term
            ansset = union(fg.nodeset[ansterm]...)
            for nset in fg.nodeset[nterm]
                if isempty(setdiff(nset, ansset))
                    push!(rmset, nset)
                end
            end
        end
        components[nterm] = setdiff(fg.nodeset[nterm], rmset)
    end
    for n in ontnodes
        rmset = Set{Set{Int}}()
        nterm = getnode(fg.graph, n).term
        for (s, e) in successors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :isa
                continue
            end
            ansterm = getnode(fg.graph, s).term
            intersect!(components[nterm], components[ansterm])
        end
        for (p, e) in predecessors(fg.graph, n)
            rel = getedge(fg.graph, e)
            if rel.relation != :isa
                continue
            end
            descterm = getnode(fg.graph, p).term
            setdiff!(components[nterm], components[descterm])
        end
    end
    return components
end
