% Create a column vector with numbers from 1 to 18
COLUMNVECTOR = (1:24)'

% Define the dimensions
x = 2;
y = 4;
z = 3;

% Reshape the column vector into a 3D array with dimensions x, y, z
reshaped_array = reshape(COLUMNVECTOR, x, y, z);

% Permute the dimensions using [2, 1, 3]
permuted_array = permute(reshaped_array, [2, 1, 3]);

% Display the original, reshaped, and permuted arrays
disp("Original Column Vector:");
disp(COLUMNVECTOR);

disp("Reshaped 3D Array:");
disp(reshaped_array);

disp("Permuted 3D Array (Dimensions swapped between x and y):");
disp(permuted_array);