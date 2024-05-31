function Triangulate(FileName)
	global N B NeighbourLeft NeighbourRight IsEar Points Triangles TriangleCount Edges EdgeCount
	
	%--------------------------------------------------
	% Reads in variables from text file
	%--------------------------------------------------
	ipf = fopen(FileName, 'r');
	N = fscanf(ipf, "%d", [1 1]);
	B = fscanf(ipf, "%d", [1 1]);
	Points = fscanf(ipf, "%f %f", [2,N]);
	fclose(ipf);

	TriangleCount = 0;				% Number of triangules currently present
	EdgeCount = 0;					% Number of Edges currently present
	Edges = zeros(2,B + 3*(N-B));	% Array of Edges made up of 2 edges, each entry corresponds to the index in the corresponding "Point" array
	Triangles = zeros(3,N-2);		% Array of triangles made up of 3 edges, each entry corresponds to the index in the corresponding "Edges" array ( +/- corresponds to orientation)

	%--------------------------------------------------
	% Ear Removal Algorithm for triangulation
	%--------------------------------------------------
	IsEar = false(1,B);
	NeighbourLeft = 0:N-1;
	NeighbourLeft(1) = B;
	NeighbourRight = 2:N+1;
	NeighbourRight(B) = 1;
	while(TriangleCount < B-3)
		for i=1:B
			if(Diagonal(NeighbourLeft(i),NeighbourRight(i))) % Check if the current vertice is an ear
				
				%Adds triangle 
				TriangleCount = TriangleCount + 1;
				Sorted = sort([NeighbourLeft(i), i, NeighbourRight(i)]);
				Triangles(1, TriangleCount) = UpdateEdges(Sorted(1),Sorted(2));
				Triangles(2, TriangleCount) = UpdateEdges(Sorted(2),Sorted(3));
				Triangles(3, TriangleCount) = UpdateEdges(Sorted(3),Sorted(1));
				
				% Removes ear by updating every points nearest neighbour
				if(TriangleCount < B-3)
					if(i>1 & i<B)
						NeighbourLeft(i+1) = NeighbourLeft(i);
						NeighbourRight(i-1) = NeighbourRight(i);
					elseif(i==1)
						NeighbourLeft(2) = NeighbourLeft(B);
						NeighbourRight(B) = NeighbourRight(2);
					else
						NeighbourLeft(1) = NeighbourLeft(B-1);
						NeighbourRight(B-1) = NeighbourRight(1);
					end
				end
				IsEar(i)=true;
			end
		end
	end
	
	% Finds the final triangle after all ears are removed
	vertices = zeros(1,3);
	CurrentSide=1;
	for i=1:B
		if(IsEar(i)==false);
			vertices(CurrentSide) = i;
			CurrentSide = CurrentSide + 1;
		end
	end

	TriangleCount = TriangleCount + 1;
	Sorted = sort(vertices);
	Triangles(1, TriangleCount) = UpdateEdges(Sorted(1), Sorted(2));
	Triangles(2, TriangleCount) = UpdateEdges(Sorted(2), Sorted(3));
	Triangles(3, TriangleCount) = UpdateEdges(Sorted(3), Sorted(1));

	figure;
	hold on;
	PlotTriangles();

	%----------------------------------------------
	% Add internal points and re-triangulates
	%----------------------------------------------

	for i = 1:N-B
		for j=1:TriangleCount
			edge1 = GetEdge(Triangles(1,j));
			edge2 = GetEdge(Triangles(2,j));
			edge3 = GetEdge(Triangles(3,j));
			
				if(LeftOn(edge1(1),  edge1(2), i+B))
					if(LeftOn(edge2(1), edge2(2), i+B)
						&& LeftOn(edge3(1),edge3(2), i+B))
							if(abs(Area2(edge1(1), edge1(2), i+B))>1.0 && abs(Area2(edge2(1), edge2(2), i+B))>1.0 % Ensures the point isn't on an edge
								&& abs(Area2(edge3(1), edge3(2), i+B))>1.0)
								UpdateTriangles(i,j);
							end
					end
				else
					if((LeftOn(edge2(1), edge2(2), i+B)==false)
						&& (LeftOn(edge3(1),edge3(2), i+B))==false)
						if(abs(Area2(edge1(1), edge1(2), i+B))>1.0 && abs(Area2(edge2(1), edge2(2), i+B))>1.0 % Ensures the point isn't on an edg
								&& abs(Area2(edge3(1), edge3(2), i+B))>1.0)
							UpdateTriangles(i,j);
						end
					end
				end
		end
	end

	figure;
	hold on;
	PlotTriangles();

	%--------------------------------------------
	% Finds optimal triangulation by edge swapping
	%--------------------------------------------
	changed=true;
	MaximumIterations = 20;
	Iteration = 0;
	while(changed)
		changed=false;
		for i=1:EdgeCount
			edge = GetEdge(i);
			if((edge(1)>=B && edge(2)>=B) || (abs(edge(2) - edge(1))>1 && ((edge(1)!=1 || edge(2)!=B) && (edge(2)!=1 || edge(1)!=B))))
				if(OptimiseTriangles(i))
					changed = true;
				end
			end
		end
		Iteration = Iteration + 1;
		if(Iteration>20)
			break
		end
	end

	figure;
	hold on;
	PlotTriangles();
	
	Triangles
	
	Edges
	
	Points
end

%---------------------------------------------------------------------------------------------------------------------------------------------

% Returns edges with order corresponding to sign of input
function e = GetEdge(i)
	global Edges
	if(i>0)
		e = Edges(:,i);
	else
		e = [Edges(2,-i), Edges(1,-i)];
	end
end

% initially determines if edge already exists if not adds edge to "Edges" array
% returns index of the edge in "Edges", sign corresponds to its orientation in "Edges" array relative to input
function e =  UpdateEdges(i,j)
	global EdgeCount Edges
	EdgeFound=false;
	for k=1:EdgeCount
		if((Edges(1,k) == i) && (Edges(2,k) == j))
				EdgeFound=true;
				e = k;
		elseif((Edges(2,k) == i) && (Edges(1,k) == j))
				EdgeFound=true;
				e = -k;
		end
	end
	if(EdgeFound==false)
		EdgeCount = EdgeCount + 1;
		Edges(1,EdgeCount) = i;
		Edges(2,EdgeCount) = j;
		e = EdgeCount;
	end
end

% Changes values of edge entry at "EdgeNumber" in "Edges" array
function ChangeEdge(EdgeNumber, i,j)
	global Edges
	Edges(1,EdgeNumber) = i;
	Edges(2,EdgeNumber) = j;
end

%Plots all current triangles in "Triangles" array
function PlotTriangles()
	global Edges N TriangleCount Points Triangles
	for  i = 1:TriangleCount
		for k=1:3
			edge = GetEdge(Triangles(k,i));
			plot([Points(1,edge(1)), Points(1,edge(2))] , [Points(2,edge(1)), Points(2,edge(2))], "b");
		end
	end
end


%---------------------------------------------------------------------------------------------------------------------------------------------
%Functions corresponding to Chapter 1 of "J. O’Rourke: Computational Geometry in C" for the "Ear Removal" algorithm
%---------------------------------------------------------------------------------------------------------------------------------------------

function a = Area2(i,j,k) % (B - A) X (C - A)
	global Points
	a = (Points(1,j) - Points(1,i))*(Points(2,k) - Points(2,i)) - (Points(1,k) - Points(1,i))*(Points(2,j) - Points(2,i));
end

function l = Left(i,j,k)
	l = Area2(i,j,k)>0;
end

function l = LeftOn(i,j,k)
	l = Area2(i,j,k)>=0;
end

function l = Collinear(i,j,k)
	l = Area2(i,j,k)==0;
end

function between=Between(i,j,k)
	if(Collinear(i,j,k)==false)
		between = false;
	else
		if(Points(1,i)!=Points(1,j))
			between =  ((Points(1,i)<=Points(1,k)) && (Points(1,k) <= Points(1,j))) || ((Points(1,i)>=Points(1,k)) && (Points(1,k) >= Points(1,j)));
		else
			between =  ((Points(2,i)<=Points(2,k)) && (Points(2,k) <= Points(2,j))) || ((Points(2,i)>=Points(2,k)) && (Points(2,k) >= Points(2,j)));
		end
	end
end

function prop=ProperIntersection(i,j,k,l)
	if(Collinear(i,j,k) ||
		Collinear(i,j,l) ||
		Collinear(k,l,i) ||
		Collinear(k,l,j))
		prop=false;
	else 
		prop = xor(Left(i,j,k), Left(i,j,l)) && xor(Left(k,l,i), Left(k,l,j));
	end
end

function inter=Intersection(i,j,k,l)
	if(ProperIntersection(i,j,k,l))
		inter=true;
	elseif(Between(i,j,k) ||
		Between(i,j,l) ||
		Between(k,l,i) ||
		Between(k,l,j))
		inter=true;
	else
		inter = false;
	end
end

function d = Diagonalie(i,j)
	global B NeighbourRight IsEar
	
	d = true;
	for k=1:B
		if((IsEar(k)==false) && (k != i) && (k != j)
			&& (NeighbourRight(k)!=i) && (NeighbourRight(k)!=j)
			&& Intersection(i,j,k, NeighbourRight(k)))
			d = false;
		end 
	end
end

function incone = InCone(i, j)
	global NeighbourLeft NeighbourRight
	if(LeftOn(i,NeighbourRight(i),NeighbourLeft(i)))
		incone = (Left(i,j,NeighbourLeft(i)) && 
				Left(j,i,NeighbourRight(i)));
	else
		incone = (LeftOn(i,j,NeighbourRight(i)) && 
				   LeftOn(j,i,NeighbourLeft(i)))==false;
	end
end

function diagonal = Diagonal(i,j)
	diagonal = InCone(i,j) && InCone(j,i) && Diagonalie(i,j);
end

%---------------------------------------------------------------------------------------------------------------------------------------------

% Adds new internal point found at index "PointNumber" of "Point" array
% Creating 3 new triangles from the new point at edges of triangle corresponding to index "TriangleNumber" of "Triangles" array
% triangle corresponding to index "TriangleNumber" is destroyed
function UpdateTriangles(InternalPointNumber, TriangleNumber)
	global B Triangles TriangleCount Edges

	TriangleCount = TriangleCount+1;
	edge = GetEdge(Triangles(2,TriangleNumber));
	Sorted = sort([edge(1),edge(2), InternalPointNumber+B]);
	Triangles(1, TriangleCount) = UpdateEdges(Sorted(1),Sorted(2));
	Triangles(2, TriangleCount) = UpdateEdges(Sorted(2),Sorted(3));
	Triangles(3, TriangleCount) = UpdateEdges(Sorted(3),Sorted(1));
		
	TriangleCount = TriangleCount+1;
	edge = GetEdge(Triangles(3,TriangleNumber));
	Sorted = sort([edge(1),edge(2), InternalPointNumber+B]);
	Triangles(1, TriangleCount) = UpdateEdges(Sorted(1),Sorted(2));
	Triangles(2, TriangleCount) = UpdateEdges(Sorted(2),Sorted(3));
	Triangles(3, TriangleCount) = UpdateEdges(Sorted(3),Sorted(1));
	
	edge = GetEdge(Triangles(1,TriangleNumber));
	Sorted = sort([edge(1),edge(2), InternalPointNumber+B]);
	Triangles(1, TriangleNumber) = UpdateEdges(Sorted(1),Sorted(2));
	Triangles(2, TriangleNumber) = UpdateEdges(Sorted(2),Sorted(3));
	Triangles(3, TriangleNumber) = UpdateEdges(Sorted(3),Sorted(1));
end



%---------------------------------------------------------------------------------------------------------------------------------------------

% Returns triangle indices of all triangles containing the edge corresponding to index "EdgeNumber" in "Edges" array 
function TriangleNumbers = FindTrianglesOfEdge(EdgeNumber)
	global Triangles TriangleCount
	
	TriangleNumbers = zeros(1,2);
	CurrentTriangle = 1;
	for i=1:TriangleCount
		for j=1:3
			if(abs(Triangles(j,i))==EdgeNumber)
				TriangleNumbers(CurrentTriangle) = i;
				CurrentTriangle= CurrentTriangle + 1;
			end
		end
	end

end
% Swaps edge corresponding to index "EdgeNumber" of "Edges" array with an edge made of opposing points
% if Opposing points correspond to a shorter length edge and is an internal edge

function changed = OptimiseTriangles(EdgeNumber)
	global Triangles Edges
	
	TriangleNumbers = FindTrianglesOfEdge(EdgeNumber);
	
	%Find Opposing to points to the "EdgeNumber" edge
	OpposingPoints = zeros(1,2);
	for i=1:2
		for k=1:3
			if(abs(Triangles(k,TriangleNumbers(i)))!=EdgeNumber)
				edge = GetEdge(Triangles(k,TriangleNumbers(i)));
				if((edge(1)!=Edges(1, EdgeNumber)) && (edge(1)!=Edges(2, EdgeNumber)))
					OpposingPoints(i) = edge(1);
				else 
					OpposingPoints(i) = edge(2);
				end
				break;
			end
		end
	end
	
	changed=false;
	if((OpposingPoints(1)!=0 && OpposingPoints(2)!=0)) % Checks if opposing points were found
		if(Length(OpposingPoints(1), OpposingPoints(2))<Length(Edges(1,EdgeNumber), Edges(2, EdgeNumber))) % Determines if edge swap is necesary for optimal triangulation
			
			temp = Edges(:,EdgeNumber);
			
			%Swap edge
			ChangeEdge(EdgeNumber, OpposingPoints(1), OpposingPoints(2));
			
			%Update triangulation
			Sorted = sort([OpposingPoints(1),OpposingPoints(2),temp(1)]);
			Triangles(1,TriangleNumbers(1)) = UpdateEdges(Sorted(1),Sorted(2));
			Triangles(2,TriangleNumbers(1)) = UpdateEdges(Sorted(2),Sorted(3));
			Triangles(3,TriangleNumbers(1)) = UpdateEdges(Sorted(3),Sorted(1));
			
			Sorted = sort([OpposingPoints(1),OpposingPoints(2),temp(2)]);
			Triangles(1,TriangleNumbers(2)) = UpdateEdges(Sorted(1),Sorted(2));
			Triangles(2,TriangleNumbers(2)) = UpdateEdges(Sorted(2),Sorted(3));
			Triangles(3,TriangleNumbers(2)) = UpdateEdges(Sorted(3),Sorted(1));
			changed=true;
		end
	end
end


% Calculates distance between points at indices, "i" and "j" from the "Point" array
function l=Length(i, j)
	global Points
	l = ((Points(1, i) - Points(1,j))^2 + (Points(2,i) - Points(2,j))^2)^(0.5);
end









