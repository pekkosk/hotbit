import sys
class Output:
    def __init__(self):
        self.txt=None
    
    def set_text(self,txt):
        """ Set the stream for text output. """
        if txt==None:
            self.txt=sys.stdout
        else:
            self.txt=open(txt,'a')
            
    def close_output(self):
        self.txt.flush()
        self.txt.close()
        
    def get_output(self):
        return self.txt
    
    def flush(self):
        self.txt.flush()