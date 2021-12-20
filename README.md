# Overview
In this final project, I have implemented the lexicographical ordering
algorithm from Ravi Sethi's 1974 paper on scheduling tasks on two cores.

I used my implementation to create graphics for my presentation,
revealing each step of the algorithm's process.

Two parts of the algorithm were not clearly described in the paper:
how one goes about implementing the graph operations, and how one
determines the level of each node.
They were both reasonably simple problems to solve, and both the
problems and solutions are described in more detail in the Auxiliary
Nodes and Determine Node Level sections.



# Ravi Sethi's paper
Although Ravi Sethi's paper includes the ordering algorithm, its
primary purpose was scheduling a known set of non-interruptable tasks
which have non-cyclic dependancies inbetween and take an equal amount
of time to run to completion; the schedule must complete all tasks
in the minimum amount of time possible.

This scheduling problem, although well solved (in linear time!),
is best suited for a dual-core embedded processor in a real-time
system, of which I have none; as such, I see little value in the
scheduling algorithm and have no motivation to implement it.


# Auxiliary Nodes
The paper operates on "auxiliary nodes" that I refer to as groups,
who are occasionally (not only at the beginning)
created with some regular (non-auxiliary) node children,
removed one-by-one no more than once,
whose children must be iterable,
and whose children must track their one parent group;
all operations must occur in constant time,
except for iterating the children and creating the node with children
(which both run in linear time with regard to the number of children).


## My Solution
The groups are in a list, where new groups are added onto the end,
and groups are never removed, so the group's indices are stable.

Each group contains a list of children node id's.

Each regular node tracks the group ID that it belongs to,
but also the index of itself in the group's children list.

To create a group with some children,
push the group to the end of the groups list,
assign the group id of each child to the index of the group in that list,
push all the children to the group's children list,
recording the index in the children list to each children's child field.

To remove a regular node from a group,
we move the node's entry in the group's children list to the end
with a swap.
The group children list elements need to be swapped,
as well as the children list index fields of both nodes.
Then, we pop the node off the end of the group's children list and clear
the group field.

The other operations should be obvious.



# Implementation
Please reference `simple.rs` for this section, as `main.rs` contains
additional code to generate my HTML+SVG+CSS slides.

The following subsections are ordered from top to bottom, as they appear
in `simple.rs`.


## Routine Label
Lines 8 through 56

This function labels each level of the graph in order of ascending level.



## Routine Partition
Lines 58 through 165


## Structs
Lines 167 through 357

`Node` is an empty struct in `simple.rs` as it's just the visual
properties of a node, which this source file isn't concerned with.

`Graph` is a directed acyclic graph, and is the input to the algorithm.
Adjacency lists are used to store both forward reverse edges.
It has a simple and moderately robust parser, implemented on lines 179
through 356.

Lines 359 through 394 define an error type used in the `Graph` parser.


## Determine Node Level
Lines 396 through 449

This function determines the level of each node in the graph in linear
time; the level is defined as one greater than the greatest child node's
level, or `0` for nodes with no children.

Ravi Sethi left the implementation details of this function as an
exercise for the reader, and so I took the liberty to implement it by
visiting all nodes with an iterative depth first search.


## fn main
Lines 451 through 512

This is the entry point to the program.
It calls parse on a hard coded DAG (feel free to modify the graph),
and then gives the parsed DAG to routine_label.
The generated labelling goes unused in `simple.rs`, but `main.rs`
visualizes it.
