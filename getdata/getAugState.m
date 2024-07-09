function augState = getAugState(state, sep)

augState = zeros(length(state(:, 1)), 13);
augState(:, 1:8) = state;

augState(sep(1,1), 9:13) = zeros(1,5);
for i = 2:sep(1,2)
    augState(i, 9:13) = state(i, 4:8) - state(i-1, 4:8);
end
augState(sep(2,1), 9:13) = zeros(1,5);
for i = sep(2,1)+1:sep(2, 2)
    augState(i, 9:13) = state(i, 4:8) - state(i-1, 4:8);
end
augState(sep(3,1), 9:13) = zeros(1,5);
for i = sep(3,1)+1:sep(3, 2)
    augState(i, 9:13) = state(i, 4:8) - state(i-1, 4:8);
end
augState(sep(4,1), 9:13) = zeros(1,5);
for i = sep(4,1)+1:sep(4, 2)
    augState(i, 9:13) = state(i, 4:8) - state(i-1, 4:8);
end
augState(sep(5,1), 9:13) = zeros(1,5);
for i = sep(5,1)+1:sep(5, 2)
    augState(i, 9:13) = state(i, 4:8) - state(i-1, 4:8);
end