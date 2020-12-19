from tkinter import *

root = Tk()
li = ['C','Python','php','html','Fortran','SQL','java']
movie = ['CSS','jQuery','Bootstrap']
listb = Listbox(root)
listb1 = Listbox(root)

for item in li:
	listb.insert(0,item)

for item in movie:
	listb1.insert(0,item)

listb.pack()
listb1.pack()
root.mainloop()
