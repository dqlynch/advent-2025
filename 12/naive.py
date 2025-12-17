shape_areas = [
    6, 5, 7, 7, 7, 7
]

with open("input.txt") as f:
    acc = 0
    for i, line in enumerate(f):
        if i < 30:
            continue
        area_str, *counts = line.strip().split(' ')
        w, h = area_str.strip(':').split('x')
        area = int(w) * int(h)
        present_area = sum([int(count) * shape_areas[i] for i, count in enumerate(counts)])

        if present_area <= area:
            acc += 1
print(acc)
