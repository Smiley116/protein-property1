import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import requests
import os
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class ProtAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("蛋白质理化性质分析工具")
        self.root.geometry("900x700")
        self.root.configure(bg='#f0f0f0')
        
        # 添加窗口置顶和焦点设置
        self.root.lift()
        self.root.attributes('-topmost', True)
        self.root.focus_force()
        self.root.update()
        
        # 设置样式
        self.setup_styles()
        self.create_widgets()
    
    def setup_styles(self):
        style = ttk.Style()
        style.configure('Title.TLabel', font=('Arial', 16, 'bold'), background='#f0f0f0')
        style.configure('Subtitle.TLabel', font=('Arial', 12, 'bold'), background='#f0f0f0')
        style.configure('Input.TLabel', font=('Arial', 10, 'bold'), background='#f0f0f0')
    
    def create_widgets(self):
        # 主容器
        main_frame = ttk.Frame(self.root, padding="15")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 标题
        title_label = ttk.Label(main_frame, text="蛋白质理化性质分析工具", style='Title.TLabel')
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # 输入方式选择
        method_frame = ttk.LabelFrame(main_frame, text="选择输入方式", padding="10")
        method_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.input_method = tk.StringVar(value="pdb")
        
        pdb_radio = ttk.Radiobutton(method_frame, text="通过 PDB ID 分析", 
                                   variable=self.input_method, value="pdb",
                                   command=self.on_method_change)
        pdb_radio.grid(row=0, column=0, sticky=tk.W, padx=(0, 20))
        
        seq_radio = ttk.Radiobutton(method_frame, text="直接输入氨基酸序列", 
                                   variable=self.input_method, value="sequence",
                                   command=self.on_method_change)
        seq_radio.grid(row=0, column=1, sticky=tk.W)
        
        # PDB输入区域
        self.pdb_frame = ttk.LabelFrame(main_frame, text="PDB ID 输入", padding="10")
        self.pdb_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        pdb_help = ttk.Label(self.pdb_frame, text="请输入4位PDB ID（例如: 1crn, 2abl, 3gb1）", 
                            background='#f0f0f0')
        pdb_help.grid(row=0, column=0, sticky=tk.W, pady=(0, 10))
        
        pdb_input_frame = ttk.Frame(self.pdb_frame)
        pdb_input_frame.grid(row=1, column=0, sticky=(tk.W, tk.E))
        
        ttk.Label(pdb_input_frame, text="PDB ID:", style='Input.TLabel').grid(row=0, column=0, sticky=tk.W)
        self.pdb_entry = ttk.Entry(pdb_input_frame, width=20, font=('Arial', 11))
        self.pdb_entry.grid(row=0, column=1, padx=(10, 0), sticky=tk.W)
        
        
        # 序列输入区域
        self.seq_frame = ttk.LabelFrame(main_frame, text="氨基酸序列输入", padding="10")
        self.seq_frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        seq_help = ttk.Label(self.seq_frame, text="请输入氨基酸序列（支持FASTA格式）", 
                            background='#f0f0f0')
        seq_help.grid(row=0, column=0, sticky=tk.W, pady=(0, 10))
        
        self.seq_text = scrolledtext.ScrolledText(self.seq_frame, height=8, width=80, 
                                                 font=('Courier', 10))
        self.seq_text.grid(row=1, column=0, sticky=(tk.W, tk.E))
        
        seq_button_frame = ttk.Frame(self.seq_frame)
        seq_button_frame.grid(row=2, column=0, pady=(10, 0))
        
        
        seq_clear_btn = ttk.Button(seq_button_frame, text="清空序列", 
                                  command=self.clear_sequence)
        seq_clear_btn.grid(row=0, column=1)
        
        # 按钮区域
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=4, column=0, columnspan=2, pady=20)
        
        self.analyze_btn = ttk.Button(button_frame, text="开始分析", 
                                     command=self.analyze_sequence,
                                     style='Accent.TButton')
        self.analyze_btn.grid(row=0, column=0, padx=(0, 10))
        
        clear_all_btn = ttk.Button(button_frame, text="清空所有", 
                                  command=self.clear_all)
        clear_all_btn.grid(row=0, column=1)
        
        # 结果显示区域
        result_frame = ttk.LabelFrame(main_frame, text="分析结果", padding="10")
        result_frame.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(10, 0))
        
        self.result_text = scrolledtext.ScrolledText(result_frame, height=20, width=80, 
                                                    font=('Courier', 9))
        self.result_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 初始显示PDB输入框
        self.on_method_change()
        
        # 配置网格权重
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(5, weight=1)
        self.pdb_frame.columnconfigure(0, weight=1)
        self.seq_frame.columnconfigure(0, weight=1)
        result_frame.columnconfigure(0, weight=1)
        result_frame.rowconfigure(0, weight=1)
    
    def on_method_change(self):
        """根据选择的输入方式显示/隐藏对应的输入框"""
        if self.input_method.get() == "pdb":
            self.pdb_frame.grid()
            self.seq_frame.grid_remove()
        else:
            self.pdb_frame.grid_remove()
            self.seq_frame.grid()
    
    def load_pdb_example(self):
        """加载PDB示例"""
        self.pdb_entry.delete(0, tk.END)
        self.pdb_entry.insert(0, "1crn")
    
    
    def clear_sequence(self):
        """清空序列输入框"""
        self.seq_text.delete(1.0, tk.END)
    
    def clear_all(self):
        """清空所有输入"""
        self.pdb_entry.delete(0, tk.END)
        self.seq_text.delete(1.0, tk.END)
        self.result_text.delete(1.0, tk.END)
    
    def analyze_sequence(self):
        """分析序列"""
        # 禁用分析按钮，防止重复点击
        self.analyze_btn.config(state='disabled')
        
        try:
            # 清空之前的结果
            self.result_text.delete(1.0, tk.END)
            
            # 显示分析中状态
            self.result_text.insert(tk.END, "分析中...\n")
            self.root.update()
            
            sequence = ""
            
            if self.input_method.get() == "pdb":
                # PDB ID 分析
                pdb_id = self.pdb_entry.get().strip()
                if not pdb_id:
                    messagebox.showerror("错误", "请输入PDB ID")
                    return
                
                sequence = self.get_sequence_from_pdb(pdb_id)
                if sequence:
                    self.result_text.insert(tk.END, f"成功获取PDB {pdb_id} 的序列\n")
                
            else:
                # 直接序列分析
                user_input = self.seq_text.get(1.0, tk.END).strip()
                if not user_input:
                    messagebox.showerror("错误", "请输入氨基酸序列")
                    return
                
                sequence = self.extract_sequence_from_input(user_input)
                if sequence:
                    self.result_text.insert(tk.END, "成功提取氨基酸序列\n")
            
            if sequence:
                # 分析序列
                self.analyze_and_display(sequence)
            else:
                self.result_text.delete(1.0, tk.END)
                self.result_text.insert(tk.END, "错误: 无法获取有效的氨基酸序列\n")
                self.result_text.insert(tk.END, "请检查输入:\n")
                if self.input_method.get() == "pdb":
                    self.result_text.insert(tk.END, "- PDB ID应为4位字符\n")
                    self.result_text.insert(tk.END, "- 请检查网络连接\n")
                else:
                    self.result_text.insert(tk.END, "- 序列应包含有效的氨基酸字符\n")
                    self.result_text.insert(tk.END, "- 序列长度应至少10个氨基酸\n")
            
        except Exception as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"分析过程中出现错误:\n{str(e)}")
        
        finally:
            # 重新启用分析按钮
            self.analyze_btn.config(state='normal')
    
    def extract_sequence_from_input(self, input_str):
        """从输入中提取序列"""
        if input_str.startswith('>'):
            lines = input_str.split('\n')
            sequence_content = ''.join(lines[1:]) if len(lines) > 1 else ""
        else:
            sequence_content = input_str
        
        return self.clean_sequence(sequence_content)
    
    def clean_sequence(self, sequence):
        """清理序列"""
        cleaned = ''.join([char for char in sequence.upper() if char.isalpha()])
        
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        sequence_chars = set(cleaned)
        
        invalid_chars = sequence_chars - valid_aa
        if invalid_chars:
            self.result_text.insert(tk.END, f"警告: 移除非标准氨基酸字符: {invalid_chars}\n")
            original_length = len(cleaned)
            cleaned = ''.join([aa for aa in cleaned if aa in valid_aa])
            self.result_text.insert(tk.END, f"已移除 {original_length - len(cleaned)} 个非标准字符\n")
        
        return cleaned
    
    def get_sequence_from_pdb(self, pdb_id):
        """从PDB获取序列"""
        if not self.is_valid_pdb_id(pdb_id):
            self.result_text.insert(tk.END, f"错误: {pdb_id} 不是有效的PDB ID\n")
            return None
        
        fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}/download"
        
        try:
            response = requests.get(fasta_url, timeout=10)
            response.raise_for_status()
            
            lines = [line.strip() for line in response.text.split('\n') if line.strip()]
            sequence_lines = [line for line in lines[1:] if "|" not in line]  
            raw_sequence = ''.join(sequence_lines).replace(' ', '').replace('\n', '')
            
            return self.clean_sequence(raw_sequence)
            
        except requests.exceptions.RequestException as e:
            self.result_text.insert(tk.END, f"下载失败: {e}\n")
            return None
    
    def is_valid_pdb_id(self, pdb_id):
        """验证PDB ID格式"""
        return len(pdb_id) == 4 and re.match(r'^[0-9a-zA-Z]{4}$', pdb_id)
    
    def analyze_and_display(self, sequence):
        """分析并显示结果"""
        prot = ProteinAnalysis(sequence)

        # 计算结果
        molecular_weight_da = prot.molecular_weight()
        molecular_weight_kda = molecular_weight_da / 1000
        isoelectric_point = prot.isoelectric_point()
        extinction_coeff_no_cys = prot.molar_extinction_coefficient()[0]
        extinction_coeff_with_cys = prot.molar_extinction_coefficient()[1]
        gravy = prot.gravy()
        aa_composition = prot.get_amino_acids_percent()

        # 显示结果
        self.result_text.delete(1.0, tk.END)
        
        # 基本属性
        self.result_text.insert(tk.END, "=" * 60 + "\n")
        self.result_text.insert(tk.END, "蛋白质理化性质分析结果\n")
        self.result_text.insert(tk.END, "=" * 60 + "\n\n")
        
        self.result_text.insert(tk.END, f"序列长度: {len(sequence)} 个氨基酸\n")
        self.result_text.insert(tk.END, f"分子量: {molecular_weight_kda:.2f} kDa ({molecular_weight_da:.2f} Da)\n")
        self.result_text.insert(tk.END, f"等电点 (pI): {isoelectric_point:.2f}\n")
        self.result_text.insert(tk.END, f"消光系数 (无二硫键): {extinction_coeff_no_cys:.0f} M⁻¹cm⁻¹\n")
        self.result_text.insert(tk.END, f"消光系数 (有二硫键): {extinction_coeff_with_cys:.0f} M⁻¹cm⁻¹\n")
        self.result_text.insert(tk.END, f"GRAVY值 (疏水性): {gravy:.3f}\n\n")
        
        # 氨基酸组成
        self.result_text.insert(tk.END, "氨基酸组成 (%):\n")
        self.result_text.insert(tk.END, "-" * 40 + "\n")
        
        aa_items = list(aa_composition.items())
        for i in range(0, len(aa_items), 4):
            line = ""
            for j in range(4):
                if i + j < len(aa_items):
                    aa, percent = aa_items[i + j]
                    line += f"{aa}: {percent*100:5.1f}%  "
            self.result_text.insert(tk.END, line + "\n")
        
        self.result_text.insert(tk.END, "\n" + "=" * 60 + "\n")
        self.result_text.insert(tk.END, "分析完成！\n")

def main():
    root = tk.Tk()
    app = ProtAnalyzerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
