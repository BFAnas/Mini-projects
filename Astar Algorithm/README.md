# A routing algorithm for finding the optimal distance between two points
To do this we constructed a graph from a map file, for this we wrote the Graph.C program
Then we implemented the A* algorithm to try to find the shortest
path (according to distance) from Basílica de Santa Maria del Mar (Plaça de Santa Maria) in
Barcelona to the Giralda (Calle Mateos Gago) in Sevilla.
We are going to test 3 variations of the A*: the firrst, with the heuristic function equal to zero,
this is the special case of the Dijkstra algorithm, the second, using the Haversine distance as the
heuristic function, and finally using the Spherical Law of Cosines.
You can find much more information in the report included: Astar.pdf
