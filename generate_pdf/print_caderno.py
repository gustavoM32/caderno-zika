import os
import re
import sys
caderno_dir = f"{sys.path[0]}/.."

out_filename = "caderno.md"
out_file = open(out_filename, 'w')

def print_code_file(path, name):
	out_file.write(f"## {name}\n")
	out_file.write("```cpp\n")
	out_file.write(open(f"{path}/{name}").read())
	out_file.write("\n```\n")
	# to break page after each file
	# out_file.write("\n<div style='page-break-after: always;'></div>\n\n")
	

def print_file(path, name, prefix=""):
	out_file.write(f"{prefix}{name}\n")
	out_file.write(open(f"{path}/{name}").read())

############################## HEADER ###################################
out_file.write("# Caderno Zikados\n\n")
out_file.write("Enrique Junchaya, ")
out_file.write("Gustavo M. Carlos, ")
out_file.write("Nathan Luiz Martins\n\n")

############################## TEMPLATE ###################################

print_file(caderno_dir, "config.md", prefix="# ")
print_file(caderno_dir, "reference.md", prefix="# ")

############################## TOPICS ###################################

out_file.write("\n")

skip_dirs = ["tools", "generate_pdf"]

def filter_topic_dirs(dir_name):
	if (dir_name[0] == '.'):
		return False
	return dir_name not in skip_dirs

topic_dirs = filter(filter_topic_dirs, list(os.walk(caderno_dir))[0][1])

def filter_files(file_name):
	if (re.search("test.cpp$", file_name)):
		return False
	return True

def beautify_topic_name(name):
	name = name.replace('_', ' ')
	name = name.capitalize()
	return name

for topic in sorted(topic_dirs):
	# page break before topic
	out_file.write("\n<div style='page-break-after: always;'></div>\n\n")	
	out_file.write(f"# {beautify_topic_name(topic)}\n")
	files = list(filter(filter_files, list(os.walk(f"{caderno_dir}/{topic}"))[0][2]))
	out_file.write(f"`{', '.join(	(files))}`\n")
	for file in files:
		print_code_file(f"{caderno_dir}/{topic}", file)
