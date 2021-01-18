include make.inc
 
TMPLIB = $(OBJDIR)/../libmytmp_$(CMODE).a
 
 
all: objs main
 
main: $(TMPLIB)
    $(LINK)   $(TMPLIB)  $(LIBDIR) $(LIBS) $(FFLAGS) -o $(EXE)
 
 
#编译所有obj文件并构造库，在make_mod后列出所有src子目录的编译命令  make -f make_sub "SUB=子目录名称"
objs: 
    make -f make_mod
    make -f make_sub "SUB=main"
    make -f make_f77 "SUB=f77"
    make -f make_clib "SUB=anis_c"
    $(ARCH) $(ARCHFLAG) $(OBJDIR)/../libmytmp_$(CMODE).a  $(OBJDIR)/*.o
 
     
#删除所有obj文件、mod文件和库，在make_mod后列出所有src子目录的编译命令  make -f make_sub "SUB=子目录名称" clean 来清理给定类别的obj文件
#这里只是举例说明可以指定clean一部分文件，在某些场合下可以采用这个方法进行较为精确的clean控制
cleanall:
    make -f make_mod  clean
    make -f make_sub "SUB=main"  clean
    make -f make_f77 "SUB=f77" clean
    make -f make_clib "SUB=anis_c" clean
    rm -f $(OBJDIR)/../libmytmp_$(CMODE).a
    rm -f $(OBJDIR)/*.mod
 
#粗暴的删除所有编译生成的文件，
clean: 
    rm -f $(OBJDIR)/../libmytmp_$(CMODE).a
    rm -f $(OBJDIR)/*.mod
    rm -f $(OBJDIR)/*.o
#创建obj文件目录
dir:
    mkdir  ../obj_$(COMPILER)
    mkdir  ../obj_$(COMPILER)/debug
    mkdir  ../obj_$(COMPILER)/release