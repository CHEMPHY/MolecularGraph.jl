#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Canvas,
    Color,
    DRAW_SETTING,
    draw2d_annot!,
    atomcolor,
    atomvisible,
    bondnotation!,
    termbondnotation!,
    ringbondnotation!,
    boundary,
    chargesign,
    atommarkup,
    atomhtml


abstract type Canvas end


struct Color
    r::Int
    g::Int
    b::Int
end


const DRAW_SETTING = Dict(
    :display_terminal_carbon => false,
    :atomcolor => Dict(
        :H => Color(0, 0, 0),
        :B => Color(128, 0, 0),
        :C => Color(0, 0, 0),
        :N => Color(0, 0, 255),
        :O => Color(255, 0, 0),
        :F => Color(0, 255, 0),
        :Si => Color(128, 64, 192),
        :P => Color(192, 0, 192),
        :S => Color(192, 192, 0),
        :Cl => Color(64, 192, 64),
        :As => Color(128, 0, 128),
        :Se => Color(128, 128, 0),
        :Br => Color(0, 192, 0),
        :I => Color(0, 128, 0)
    ),
    :defaultatomcolor => Color(0, 192, 192)
)


function draw2d_annot!(mol::VectorMol,
                       setting=copy(DRAW_SETTING); recalculate=false)
    if haskey(mol, :Coords2D) && !recalculate
       return
    end
    topology!(mol, recalculate=recalculate)
    elemental!(mol, recalculate=recalculate)
    if nodetype(mol) === SDFileAtom
        mol[:Coords2D] = zeros(Float64, atomcount(mol), 2)
        for (i, a) in nodesiter(mol)
            mol[:Coords2D][i, :] = a.coords[1:2]
        end
    else
        mol[:Coords2D] = coords2d(mol)
    end
    mol[:AtomColor] = atomcolor(setting).(mol[:Symbol])
    mol[:AtomVisible] = atomvisible(
        setting[:display_terminal_carbon]).(mol[:Symbol], mol[:Degree])
    bondnotation!(mol)
    return
end


function atomcolor(setting)
    symbol -> get(setting[:atomcolor], symbol, setting[:defaultatomcolor])
end

function atomvisible(termc)
    (symbol, numb) -> numb == 0 || symbol != :C || (termc && numb == 1)
end


function bondnotation!(mol::VectorMol)
    if nodetype(mol) === SDFileAtom
        mol[:BondNotation] = getproperty.(edgevector(mol), :notation)
    else
        mol[:BondNotation] = zeros(Int, bondcount(mol))
    end
    termbondnotation!(mol)
    ringbondnotation!(mol)
    return
end


function termbondnotation!(mol::VectorMol)
    for a in findall(mol[:Degree] .== 1)
        termbond = pop!(neighboredgekeys(mol, a))
        if mol[:BondOrder][termbond] == 2
            mol[:BondNotation][termbond] = 2
        end
    end
    return
end


function ringbondnotation!(mol::VectorMol)
    rings = mol.annotation[:Topology].rings
    coords = mol[:Coords2D]
    for ring in sort(rings, by=length, rev=true)
        vtcs = [vec2d(coords[n, :]) for n in ring]
        cw = isclockwise(vtcs)
        if cw === nothing
            continue
        end
        succ = Dict()
        ordered = cw ? ring : reverse(ring)
        for n in 1:length(ring)-1
            succ[ordered[n]] = ordered[n+1]
        end
        succ[ordered[end]] = ordered[1]
        rsub = nodesubgraph(mol.graph, ring)
        for (i, e) in edgesiter(rsub)
            if mol[:BondOrder][i] != 2
                continue
            end
            if succ[e.u] == e.v
                mol[:BondNotation][i] = 1
            end
        end
    end
    return
end


function boundary(mol::VectorMol, coords)
    (left, right) = extrema(coords[:, 1])
    (bottom, top) = extrema(coords[:, 2])
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for bond in edgevector(mol)
        u = vec2d(coords[bond.u, :])
        v = vec2d(coords[bond.v, :])
        d = norm(v - u)
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        if long > 0.0001
            unit = long / sqrt(atomcount(mol))
        else
            unit = 1
        end
    else
        unit = median(dists) # Median bond length
    end

    return (top, left, width, height, unit)
end


function chargesign(charge)
    if charge == 0
        return ""
    end
    sign = charge > 0 ? "+" : "–" # en dash, not hyphen-minus
    num = abs(charge)
    return num > 1 ? string(num, sign) : sign
end


function atommarkup(symbol, charge, hcount, direction,
                    substart, subend, supstart, supend)
    if hcount == 1
        hs = "H"
    elseif hcount > 1
        hs = string("H", substart, hcount, subend)
    else
        hs = ""
    end
    chg = charge == 0 ? "" : string(supstart, chargesign(charge), supend)
    txt = direction == :left ? [chg, hs, symbol] : [symbol, hs, chg]
    return string(txt...)
end


function atomhtml(symbol, charge, hcount, direction)
    atommarkup(symbol, charge, hcount, direction,
               "<sub>", "</sub>", "<sup>", "</sup>")
end
