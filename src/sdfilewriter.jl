#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export



import Base: parse, iterate, eltype, IteratorSize, IteratorEltype


struct SDFileReader
    lines::Base.EachLine
    parser::Function
end


"""
    sdfilereader(file::IO)

Read SDFile data from input stream and return a lazy iterator which
yields molecule objects.

`sdfilereader` does not stop and raise errors when an erroneous or incompatible
SDFile block is read but produces an error message and yields an empty molecule.
If this behavior is not desirable, you can use the customized supplier function
instead of default supplier `nohaltsupplier`

"""
sdfilewriter(mols, file::IO) = sdfilewriter(mols)


"""
    sdfilereader(path::AbstractString)

Read a SDFile and return a lazy iterator which yields molecule objects.
"""
sdfilewriter(mols, path::AbstractString) = sdfilewriter(open(path))


function molblock(mol::MolGraph)
    # TODO: convert to vectormol
    # TODO: 2D coords for SmilesMol

    lines = String[
        "",
        "MolecularGraph.jl version $(version())",
        ""
    ]
    chiral_flag = 0 # TODO: deplicated?
    acnt = @sprintf "%3d" nodecount(mol)
    bcnt = @sprintf "%3d" edgecount(mol)
    cflag = @sprintf "%3d" edgecount(chiral_flag)
    push!(lines, "$(acnt)$(bcnt)  0  0$(cflag)  0  0  0  0  0999 V2000")

    (atoms, idx_table) = atomblock(mol)
    append!(lines, atoms)
    bonds = bondblock(mol, idx_table)
    if bonds
        append!(lines, bonds)
    end
    props = propblock(mol, idx_table)
    if props:
        lines.extend(props)
    end
    push!(lines, "M  END")
    if sdfile
        append!(datablock(mol))
    end
    push!(lines, "")
    return join(lines, "\n")
end


function atomblock(mol::VectorMol)
    lines = String[]
    for (i, atom) in nodesiter(mol)
        (_x, _y, _z) = atom.coords
        x = @sprintf "%010.4f" _x
        y = @sprintf "%010.4f" _y
        z = @sprintf "%010.4f" _z
        sym = @sprintf "%-3d" string(atom.symbol)
        push!(lines, "$(x)$(y)$(z) $(sym) 0  0  0  0  0  0  0  0  0  0  0  0")
    end
    return lines
end


function bondblock(mol::VectorMol)
    stereo_conv = {
        1: {0: 0, 1: 1, 3: 4, 2: 6},
        2: {0: 0, 1: 1, 3: 3, 2: 6}
    }
    lines = []
    for (i, bond) in edgesiter(mol)
        u = @sprintf "%3d" bond.u
        v = @sprintf "%3d" bond.v
        order = @sprintf "%3d" bond.order
        bond_line.append("{:>3}".format(idx_table[u]))
        bond_line.append("{:>3}".format(idx_table[v]))
        bond_line.append("{:>3}".format(b.order))
        # TODO: SmilesBond
        bond_line.append("{:>3}".format(b.notation))
        push!(lines, "$(u)$(v)$(order)$(stereo)  0  0  0")
    return lines
end
