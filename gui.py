from tkinter import *
from tkinter import filedialog
from sanger import *

# The main window
root = Tk()
root.geometry("1000x1000")

# Title
title = Label(root, text="Sanger Sequencing")
title.pack()

label_file_explorer = Label(root, 
                            text = "Choose File or Folder",
                            width = 100, height = 4, 
                            fg = "blue")
label_file_explorer.pack()

# Second Prompt for Files
def filePrompt():
    prompt = Toplevel(root)

    title1 = Label(prompt, text = "What are you trying to run?", font='Helvetica 18 bold')

    file_button = Button(prompt, text = "File (.ab1, .zip)", 
                         command=lambda: [browseFiles(), prompt.destroy()], 
                         height=5, width=30, bg="yellow")
    dir_button = Button(prompt, text = "Directory (Folders)", 
                        command=lambda: [browseFolders(), prompt.destroy()], 
                        height=5, width = 30, bg="yellow")
    
    text1 = Label(prompt, text="Files are for single .ab1 files or zipped directories")
    text2 = Label(prompt, text="Directories are only for folders with .ab1 files (will not unzip)")

    title1.pack()
    text1.pack()
    file_button.pack()
    text2.pack()
    dir_button.pack()

# Buttons
b1 = Button(root, text="Choose", width=50, command=filePrompt)
b1.pack()

# Opens up File explorer and asks to choose a folder or file
selected = ""
def browseFiles():
    global selected

    file_path = filedialog.askopenfilename(filetypes=(("Zip files",'.zip'), ("Ab1 Files", '.ab1')))
    if file_path:
        selected = os.path.abspath(file_path)
    
    label_file_explorer.configure(text="File Opened: " + selected, fg="blue")

def browseFolders():
    global selected
    path = filedialog.askdirectory()
    if path:
        selected = os.path.abspath(path)

    label_file_explorer.configure(text="Folder Opened: " + selected, fg="blue")

# Runs Sanger on the selection
def checkFolder():
    global selected

    good_seqs = []
    directory = selected
    base = os.path.basename(selected)

    if not selected:
        label_file_explorer.configure(text="Nothing Selected!", fg="red")
        return
    
    else:
        if base.endswith(".zip"):
            base = base[:-4]
            directory = os.path.dirname(selected)
            extracted = extractZip(selected)
            files = findAb1Files(extracted)

        elif selected.endswith('.ab1'):
            base = base[:-4]
            directory = os.path.dirname(selected)
            files = [selected]

        else:
            files = findAb1Files(selected)

        print("Selected: " + selected)
        fastq_folder = convertAb1Fastq(files, (directory, base))

        if fastq_folder == "empty":
            return

        for fastq in os.listdir(fastq_folder):
            header, seq, qual, avg_qual, count = \
                calcFastqQuality(os.path.join(os.path.abspath(fastq_folder), fastq))
            
            print(f"{header}\nAverage Quality: {avg_qual:.2f}\n# of uncertain: {count}")

            if avg_qual > 20:
                good_seqs.append(header + ": " + str(avg_qual))

    print("Single Strains:")
    print(good_seqs)

b2 = Button(root, text="Run", command=checkFolder)
b2.pack()

root.mainloop()