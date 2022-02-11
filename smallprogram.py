from getseq import *
import data_analysis as dt
import dp as dp
from tkinter import *

class Smallprogram:
    def __init__(self, model):
        self.model = model
        self.state = 'off'

        self.root = Tk()
        self.root.title(self.model)
        self.lb_state = Label(self.root, text=self.state, fg='blue', bg='yellow')
        self.lb_state.grid()
        self.bn_state = Button(self.root, text='On/Off', command=self.on_off)
        self.bn_state.grid()

        self.bn_fun1 = Button(self.root, text='Sequence Alignment ', command=self.fun1)
        self.bn_fun1.grid(row=0, column=1, columnspan=2, rowspan=2)
        self.bn_fun2 = Button(self.root, text='Analyze the law of data', command=self.fun2)
        self.bn_fun2.grid(row=0, column=3, columnspan=2, rowspan=2)
        self.bn_fun3 = Button(self.root, text='Generating random sequences', command=self.fun3)
        self.bn_fun3.grid(row=0, column=5, columnspan=2, rowspan=2)

        self.lb_seq1 = Label(self.root, text="First sequence")
        self.lb_seq1.grid(row=3, column=1)
        self.Et_seq1 = Entry(self.root, bd=2)
        self.Et_seq1.grid(row=4, column=1)
        self.lb_seq2 = Label(self.root, text="Second sequence")
        self.lb_seq2.grid(row=5, column=1)
        self.Et_seq2 = Entry(self.root, bd=2)
        self.Et_seq2.grid(row=6, column=1)


        self.lb_n = Label(self.root, text="Number of sequences")
        self.lb_n.grid(row=3,column=5)
        self.Et_n = Entry(self.root, bd=2)
        self.Et_n.grid(row=4,column=5)
        self.lb_m = Label(self.root, text="Sequence length")
        self.lb_m.grid(row=5,column=5)
        self.Et_m = Entry(self.root, bd=2)
        self.Et_m.grid(row=6,column=5)

        self.notion1 = Label(self.root, text="Welcome to use this tool for sequence alignment\nCheck the source code for help", fg='red')
        self.notion1.grid(row=0, column=7, rowspan=2)
        self.notion2 = Label(self.root, text="No response to the program")
        self.notion2.grid(row=2, column=7, rowspan=5)

        self.root.mainloop()

    def on_off(self):
        if self.state == 'off':
            self.state = 'on'
        else:
            self.state = 'off'
        self.lb_state['text'] = self.state

    def fun1(self):
        """实现第一个功能"""
        if self.state == 'on':
            seq1=self.Et_seq1.get()
            seq2=self.Et_seq2.get()
            (max_score,s1,s2)=dp.main(seq1, seq2)
            self.notion2['text']="The max score is：" + str(max_score) + '\n' + "The best matching method is：" + '\n' + s1 + '\n' + s2 + '\n' + "Sequence alignment completed！"

    def fun2(self):
        """实现第二个功能"""
        if self.state == 'on':
            outcome=dt.main()
            self.notion2['text'] = "KS results are as follows：\n"+str(outcome)+"Draw and save the img\nData analysis completed！\nYou can see the outcome file for details"

    def fun3(self):
        """实现第三个功能"""
        if self.state == 'on':
            n=self.Et_n.get()
            m=self.Et_m.get()
            n=int(n)
            m=int(m)
            seqs=[]
            for i in range(n):
                seqs.append(getseq(m))

            SEQ = open("seqs.txt", 'w')
            for i in seqs:
                SEQ.write(''.join(i))
                SEQ.write('\n')

            self.notion2['text'] = "The first sequence produced is：\n" + ''.join(seqs[0]) + "\nOther sequences have been saved in the file\nGenerating random sequences completed！"


smallprogram = Smallprogram('Gene comparison')