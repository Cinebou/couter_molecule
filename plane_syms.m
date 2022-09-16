
syms ax bx by cx cy cz x y z real
a = [ax 0 0];
b = [bx by 0];
c = [cx cy cz];

XYZ = [x y z];
x_min_plane = dot(cross(b,c),XYZ); 
x_max_plane = dot(cross(b,c),(XYZ-a)); 

y_min_plane = dot(cross(a,c),XYZ);
y_max_plane = dot(cross(a,c),(XYZ-b));