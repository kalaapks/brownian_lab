from brownian_lab import donsker_walk
#it is coded inside __init__ in this manner

t, W = donsker_walk(20, seed=43)
# print("Time grid t:\n", t)
# print("\nWalk values W:\n", W)
print("\n(Time, Walk) coordinates:\n")
for time_pt, walk_val in zip(t, W):
    print(f"({time_pt:.2f}, {walk_val:.4f})")