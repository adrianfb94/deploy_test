import os

def find_files(filename, search_path):
   result = []

# Wlaking top-down from the root
   for root, dir, files in os.walk(search_path):
      if filename in files:
         result.append(os.path.join(root, filename))
   return result


result = find_files(".bashrc","/home")

# find_log = open('find.log', 'w')

for r in result:
   print(r)

#   find_log.write(r)
#find_log.close()


