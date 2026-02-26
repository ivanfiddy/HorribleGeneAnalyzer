import time
from io import TextIOWrapper
import pprint
import matplotlib.pyplot as plt
import numpy
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox

start = time.time()


class GeneCode:
    def __init__(self, context):
        self.PTG_hash = {
            'phe': ['uuu', 'uuc'],
            'leu': ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
            'ile': ['auu', 'auc', 'aua'],
            'met': ['aug'],
            'val': ['guu', 'guc', 'gua', 'gug'],
            'ser': ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
            'pro': ['ccu', 'ccc', 'cca', 'ccg'],
            'thr': ['acu', 'acc', 'aca', 'acg'],
            'ala': ['gcu', 'gcc', 'gca', 'gcg'],
            'tur': ['uau', 'uac'],
            'stop': ['uaa', 'uag', 'uga'],
            'his': ['cau', 'cac'],
            'gln': ['caa', 'cag'],
            'asn': ['aau', 'aac'],
            'lys': ['aaa', 'aag'],
            'asp': ['gau', 'gac'],
            'glu': ['gaa', 'gag'],
            'cys': ['ugu', 'ugc'],
            'trp': ['ugg'],
            'arg': ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
            'gly': ['ggu', 'ggc', 'gga', 'ggg']
        }

        self.CTA_hash = {}
        self.codons = []
        self.og_gene = None

        self.__context(context)
        self.__ctp_convertor()
        self.__codons()
        pass

    def __context(self, context, typ=TextIOWrapper):

        print('Beginning file reading')
        if isinstance(context, str) and typ == str:
            self.og_gene = numpy.array(list(context.lower().replace('\n', '')))
        elif isinstance(context, TextIOWrapper) and typ == TextIOWrapper:
            cont = context.read()
            if cont[0] == '>':
                self.og_gene = numpy.array(list(cont.lower()[cont.index('\n'):].replace('\n', '')))
            else:
                self.og_gene = numpy.array(list(cont.lower().replace('\n', '')))
        elif isinstance(context, list) and typ == list:
            self.og_gene = numpy.array([i.lower() for i in context])
        else:
            raise ValueError("The argument entered was neither a string or a file that can be read")
        print('File reading complete, self.og_gene saved')

    def __ctp_convertor(self):
        for i, j in self.PTG_hash.items():
            for x in j: self.CTA_hash[x] = i
            pass
        pass

    def __codons(self):
        pre_codons = self.og_gene.copy()
        translation_dict = {'a': 'u',
                            't': 'a',
                            'g': 'c',
                            'c': 'g'}
        for i in range(len(pre_codons)):
            try:
                pre_codons[i] = translation_dict[pre_codons[i]]
            except Exception as e:
                print(i)
                print(pre_codons[i])
                print('Error in __codons\n' + str(e))
        print('Saving self.codons')
        self.codons = numpy.array([''.join(pre_codons[i:i + 3]) for i in range(0, len(pre_codons), 3)])
        print('Self.codons complete')

    def BaseRatio(self):
        g_c = numpy.count_nonzero(self.og_gene == 'g') + numpy.count_nonzero(self.og_gene == 'c')
        t_a = numpy.count_nonzero(self.og_gene == 't') + numpy.count_nonzero(self.og_gene == 'a')
        return f'Ratio GC/TA = {g_c / t_a}'

    def Translate(self):
        seq = []
        for i in self.codons:
            try:
                seq.append(self.CTA_hash[i])
            except Exception as e:
                print('Error in def translate, at element' + str(i) + " and exception is: " + str(e))
                seq.append('___')
        return seq
        pass

    def CodonCounter(self):

        counter = {i: 0 for i in self.CTA_hash.keys()}

        for i in self.codons:
            try:
                counter[i] += 1
            except:
                pass

        return pprint.pformat(counter, width=100)

    def ReverseComplement(self):
        rcomplement = []
        transdict = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
        for i in self.og_gene[::-1]:
            rcomplement += [transdict[i]]
        return rcomplement


class AminoAcids(GeneCode):

    def __init__(self, inheritance):
        super().__init__(inheritance)
        pass

    def Potential(self, start=0, stop=0, limit=None, overlap=True):

        if not stop: stop = len(self.codons) - 1
        amino_sequences = []
        codon_index = []

        def __indextoamino(index):
            try:
                return self.CTA_hash.get(self.codons[index])
            except Exception as e:
                print('ERROR WITH CTP:\n' + str(e) + f"\nError in index: {index}, giving codon: {self.codons[index]}")
                return '___'

        while start <= stop:

            if __indextoamino(start) == 'met':
                xindex = [start]
                seq = [['met']]
                start += 1

                while start <= stop and __indextoamino(start) != 'stop':
                    if overlap and __indextoamino(start) == 'met': seq.append([]); xindex.append(start)
                    for j in seq: j.append(__indextoamino(start))
                    start += 1
                    pass

                if len(xindex) > 1:
                    for j in range(len(xindex)):
                        amino_sequences.append(seq[j])
                        codon_index.append((xindex[j], start))
                else:
                    amino_sequences.append(seq[0])
                    codon_index.append((xindex[0], start))

                if limit and len(amino_sequences) >= limit:
                    return amino_sequences, codon_index
                else:
                    start += 1; continue

            start += 1

        return amino_sequences, codon_index

    def ASCII_PotentalGenes(self, start=0, stop=0, dimensions='170*25'):

        if stop == 0: stop = len(self.codons) - 1

        dim = [int(i) for i in (dimensions.split('*'))]
        chart = list('-' * (dim[0] * dim[1]))
        indexes = self.Potential(start, stop)[1]
        frac = (indexes[-1][-1] - indexes[0][0]) / (dim[0] * dim[1])

        for i in indexes:
            running_pos = round((i[1] - i[0]) / frac)
            start_pos = round(i[0] / frac)
            chart[start_pos:start_pos + running_pos] = ['*'] * running_pos

        chart.insert(0, '1: ')
        for newl in range(dim[1]): chart.insert((newl + 1) * dim[0], f'\n{newl + 2}: ')

        print(f"""
        Visualization of potential coding regions in the genetic code
        Each '*' stands for approximately {round(frac, 2)} coding triplet
        Graph dimensions = {dimensions}, equal to {dim[0] * dim[1]} characters
        From / to CODONS indexes = {start} to {stop} out of total {len(self.codons)} codons ({round((stop - start) / len(self.codons) * 100, 2)}%)


{''.join(chart)}""")
        return

    def BandsGraph(self, start=0, stop=0, y_indexes=10, limit=None, overlap=False, toplevelroot=None):
        if not stop: stop = len(self.codons) - 1

        og_sequences = self.Potential(start, stop, limit, overlap)
        indexes = og_sequences[1]
        per_row = round(len(indexes) / y_indexes)
        row_format = [indexes[i * per_row:per_row * (i + 1)] for i in range(y_indexes)]
        data = [[]]

        for key, value in enumerate(row_format):
            diff = value[0][0]
            row_length = value[-1][-1] - value[0][0]
            for j in value:
                n1 = float(((j[0] - diff) / (row_length)) * 10)
                n2 = float(((j[1] - j[0]) / (row_length)) * 10)
                data[-1].append((n1, n2))
            data.append([])

        del data[-1]

        fig, ax = plt.subplots()
        y_pos = -0.2
        for key, i in enumerate(data):
            ax.broken_barh(i, (y_pos, 0.8), picker=True).set_gid(row_format[key]);
            y_pos += 1
        ax.set_yticks(range(y_indexes))
        ax.invert_yaxis()
        ax.set_xlim(0, 10)

        def raw_click(event):
            position = event.artist.get_gid()[event.ind[0]]
            amino_acid = og_sequences[0][og_sequences[1].index(position)]
            print(amino_acid)
            print(position)
            pass

        def gui_click(event):
            position = event.artist.get_gid()[event.ind[0]]
            amino_acid = '-'.join(og_sequences[0][og_sequences[1].index(position)])
            result_txt.delete('1.0', 'end')
            result_txt.insert('1.0', f"Amino Acid: \n{amino_acid}\n\nIndex: \n{position}")
            bandwindow.update_idletasks()

        if isinstance(toplevelroot, object):
            bandwindow = Toplevel(toplevelroot)
            bandwindow.title("Plot information")
            bandwindow.geometry('500x200')
            bandframe = ttk.Frame(bandwindow, relief='groove')
            result_txt = Text(bandframe)
            result_txt.insert('1.0', 'Press on a band to receive the ORF\'s information')

            bandwindow.columnconfigure(0, weight=1)
            bandwindow.rowconfigure(0, weight=1)
            bandframe.columnconfigure(0, weight=1)
            bandframe.rowconfigure(0, weight=1)
            bandframe.grid(row=0, column=0, sticky=(N, W, E, S))
            result_txt.grid(row=0, column=0, sticky=(N, W, E, S), padx=5, pady=5)
            fig.canvas.mpl_connect('pick_event', gui_click)
            bandwindow.after(ms=1000, func=plt.show)
            bandwindow.mainloop()
            pass
        else:
            fig.canvas.mpl_connect('pick_event', raw_click)
            plt.show()

        pass

    pass


class GUI():

    def __init__(self):
        self.parent = None
        self.selected_file = None
        self.root = self.__initilizer()
        self.root.mainloop()

    def __initilizer(self):

        # Tasks; Add theme/style, use numpy for better workings, convert proteins to one characters

        # Root GUI
        root = Tk()
        root.title("Human Gene Analysis")
        root.iconbitmap(True, "DNA icon.ico")
        root.geometry('500x200')

        question_mark = PhotoImage(file='questionmark.png').subsample(16, 16)
        folder = PhotoImage(file='folder.png').subsample(20, 20)

        dnadict = {'Original gene sequence': self.__GTC, 'Complementary reverse sequence': self.__comprevseq,
                   'Ratio of bases': self.__baseratio, 'Triple codon sequence': self.__triplecodseq,
                   'Frequency of codons': self.__codonfreq}
        genedict = {'Translate sequence': self.__aminotrans, 'All ORFs': self.__allorfs,
                    'Graph bands ORFs': self.__bandorfs}

        # Frames and widgets
        back_frame = ttk.Frame(root, padding=[10] * 4, relief='groove')
        footnote = ttk.Frame(root, padding=[2] * 4, relief='solid', height=30)
        author_label = ttk.Label(footnote, text='Author: Yan Fidelskiy', font=('Default', 8), anchor='n')
        help_button = ttk.Button(back_frame, text='Info', image=question_mark, compound='right',
                                 command=self.__infobutton);
        help_button.image = question_mark
        select_file = ttk.Button(back_frame, text='Select Fasta file', image=folder, compound='right',
                                 command=lambda: self.__setfilename(label_selfile));
        select_file.image = folder
        dna_analysis = ttk.Combobox(back_frame, textvariable=StringVar(), state='readonly', width=30)
        dna_analysis['values'] = list(dnadict.keys())
        dna_analysis.set("Select analysis method")
        gene_analysis = ttk.Combobox(back_frame, textvariable=StringVar(), state='readonly', width=30)
        gene_analysis['values'] = list(genedict.keys())
        gene_analysis.set("Selected analysis method")
        dna_label = ttk.Label(back_frame, text='DNA analysis:')
        gene_label = ttk.Label(back_frame, text='Gene analysis:')
        label_selfile = ttk.Label(back_frame, text='No file selected', borderwidth=1, relief='solid', justify=CENTER,
                                  padding=[5] * 4)

        # Gridding
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)

        # - idk if this is working
        back_frame.columnconfigure(1, weight=1)
        back_frame.rowconfigure(1, weight=1)
        footnote.columnconfigure(1, weight=1)
        # back_frame.grid_rowconfigure(0, weight=1)
        # back_frame.grid_columnconfigure(0, weight=1)

        # back_frame.columnconfigure(1, weight=1)
        # back_frame.rowconfigure(1, weight=1)
        # footnote.columnconfigure(1, weight=1)

        back_frame.grid(column=0, row=0, sticky=(N, W, E, S), padx=2, pady=2)
        footnote.grid(row=2, sticky=(N, W, E, S))
        author_label.grid(column=3, sticky='e')
        help_button.grid(row=0, column=1, sticky='s')
        select_file.grid(row=1, column=1, sticky='n')
        dna_analysis.grid(row=1, column=0, sticky='n')
        gene_analysis.grid(row=3, column=0, sticky='n')
        dna_label.grid(row=0, column=0, sticky='s')
        gene_label.grid(row=2, column=0, sticky='s')
        label_selfile.grid(row=2, column=1)

        # Events
        dna_analysis.bind("<<ComboboxSelected>>", lambda x: dnadict[dna_analysis.get()](
            root) if self.selected_file != None else self.__nofile())
        gene_analysis.bind("<<ComboboxSelected>>", lambda x: genedict[gene_analysis.get()](
            root) if self.selected_file != None else self.__nofile())

        return root
        pass

    def __infobutton(self):

        infoframe = ttk.Frame(self.root, relief='groove')
        infoframe.columnconfigure(1, weight=1)
        infoframe.rowconfigure(1, weight=1)
        infoframe.grid(row=0, column=0, padx=2, pady=2, sticky=(N, W, E, S))

        helplabel = ttk.Label(infoframe, font=('Default', 15), text='*finish this explaining how program works*',
                              anchor='center')
        helplabel.grid(row=1, column=1)

        leavebutton = ttk.Button(infoframe, text='X', command=infoframe.destroy)
        leavebutton.tkraise()
        leavebutton.grid(row=1, column=1, sticky='ne')

        pass

    def __setfilename(self, label):
        self.selected_file = filedialog.askopenfilename().replace('/', '\\')
        label['text'] = ("File selected:\n" + self.selected_file.split('\\')[-1])

        # Proress bar setup

        # Loading file

        with open(self.selected_file, 'r') as f:
            self.parent = AminoAcids(f)
            pass

        ### Loads fast enough
        # loading = Toplevel(self.root, relief='groove')
        # loading.rowconfigure(1, weight=1)
        # loading.columnconfigure(1, weight=1)
        # loading_bar = ttk.Progressbar(loading, orient='horizontal', length=100, mode='indeterminate')
        #
        # loading_bar.grid(row=1, column=1)
        #
        # def loading_proc():
        #     loading_bar.start(1000)
        #     with open(self.selected_file, 'r') as f: self.parent = AminoAcids(f)
        #     loading.destroy()
        #     pass
        #
        # loading.after(1, loading_proc)
        # loading.mainloop()

        pass

    def __nofile(self):
        messagebox.showerror(title='Error', message='Select Fasta file first')
        pass

    def __GTC(self, root):

        self.__displaylist(self.parent.genelist)

        pass

    def __comprevseq(self, root):

        self.__displaylist(self.parent.ReverseComplement())

        pass

    def __baseratio(self, root):

        self.__displaylist(self.parent.BaseRatio())

        pass

    def __triplecodseq(self, root):

        self.__displaylist(self.parent.codons, sep='-', chars=3)

        pass

    def __codonfreq(self, root):

        self.__displaylist(list(self.parent.CodonCounter()))

        pass

    def __aminotrans(self, root):

        self.__displaylist(self.parent.Translate(), sep='-', chars=3)

        pass

    def __allorfs(self, root):
        seq = []
        aminos, pos = self.parent.Potential()
        for i in range(len(pos)):
            seq.append(f"""Amino Acid:\n{'-'.join(aminos[i])}\n\nPosition: {str(pos[i])}""")

        self.__displaylist(seq, tpf=1)
        pass

    def __bandorfs(self, root):

        self.parent.BandsGraph(toplevelroot=self.root, overlap=False)

        pass

    def __displaylist(self, result_list, chars=1, tpf=500, sep=''):

        display = Toplevel(self.root)
        displayframe = ttk.Frame(display, relief='groove', height=400, width=700)
        result_txt = Text(displayframe)

        display.columnconfigure(0, weight=1)
        display.rowconfigure(0, weight=1)
        displayframe.columnconfigure(0, weight=1)
        displayframe.rowconfigure(0, weight=1)

        if tpf > len(result_list):
            self.tpf = len(result_list) - 1
        else:
            if chars > 1:
                self.tpf = round(tpf / chars)
            else:
                self.tpf = tpf

        def info_length():
            return numpy.floor(len(result_list) / self.tpf)

        result_txt.insert('1.0', f'Position: 0, 0%\n\n{sep.join(result_list[:self.tpf])}')

        result_txt['state'] = 'disabled'

        def update_displaytxt(val):
            current_pos.set(val)
            pos = round(int(val.split('.')[0]) * self.tpf)
            result_txt['state'] = 'normal'
            result_txt.delete('1.0', 'end')
            result_txt.insert('1.0',
                              f'Position: {pos}, {round((pos / len(result_list) * 100), 2)}%\n\n{sep.join(result_list[pos:pos + self.tpf])}')
            result_txt['state'] = 'disabled'
            display.update_idletasks()
            pass

        def update_tpf():
            try:
                self.tpf = int(string_tpf.get())
            except Exception as e:
                messagebox.showerror("Invalid input", "Enter valid integer")
                return
            update_displaytxt(current_pos.get())
            scale.config(to=info_length())
            display.update_idletasks()
            pass

        label_tpf = ttk.Label(displayframe, text="Amount of elements displayed per frame: ")
        string_tpf = StringVar()
        current_pos = StringVar();
        current_pos.set('0')
        button_tpf = ttk.Button(displayframe, command=update_tpf, text='Apply')
        entry_tpf = ttk.Entry(displayframe, textvariable=string_tpf)
        scale = ttk.Scale(displayframe, orient='vertical', length=400, from_=0, to=info_length(),
                          command=update_displaytxt)

        scale.grid(row=0, column=1)
        displayframe.grid(row=0, column=0, sticky=(N, W, E, S))
        result_txt.grid(row=0, column=0, sticky=(N, W, E, S), padx=25, pady=10)
        label_tpf.grid(row=1, column=0, sticky=(N, W, E, S))
        entry_tpf.grid(row=1, column=1, sticky=(N, W, E, S))
        button_tpf.grid(row=1, column=2, sticky=(N, W, E, S))

        display.mainloop()
        del self.tpf
        pass

    pass


### Executables here 

GUI()

### File information
# chromosome 21.fasta is an example fasta file that can be analyed with this program

print(f'Code took: {time.time() - start} seconds to finish')

# Exit code 0 meaning successful program execution
exit(0)
