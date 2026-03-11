# Elongate the provided image horizontally and save the result
from PIL import Image

input_path = "./h3_globe.png"
output_path = "./elongated_h3_globe.png"

img = Image.open(input_path)
w, h = img.size

# Stretch horizontally by 1.5x while keeping height the same
new_w = int(w * 1.3)
elongated = img.resize((new_w, h), Image.BICUBIC)

elongated.save(output_path)

print("Original size:", (w, h))
print("New size:", (new_w, h))
print("Saved to:", output_path)