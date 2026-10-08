真实路径是 /Work_bio/dropbox/Dropbox_Data/Wei_Qi，和 files/ 在同一块盘 /dev/sdb1，bind mount 可以直接用，不占额外空间。

需要你执行的 sudo 命令（我没有 sudo 权限，这三行是整套操作）：


# 第一行是用 root 权限新建一个空目录，同时指定它的属主、属组和权限。

# 部分	含义
# sudo	用 root 权限执行。files/ 属于 client_sftp，gao 在里面建不了子目录，所以需要它。
# install -d	install 的目录模式，相当于 mkdir -p 加上设置属主和权限，一步完成。
# -o gao	新目录的属主是 gao。这样以后你不用 sudo 也能管理这个目录。
# -g sftp_users	属组是 sftp_users，和 files/ 的属组一致。
# -m 755	权限为 rwxr-xr-x：gao 可写，其他人（包括 client_sftp）只能读和进入。
# 最后的路径	要建的目录，就是 bind mount 的挂载点。
# 它只建一个空目录，不改 files/ 本身的权限，也不碰任何已有文件。

# 这个目录只是"挂载点"。bind mount 挂上以后，里面显示的是 Dropbox 里 Wei_Qi 的内容，数据一个字节都没复制。

# 你刚打开的 /home/gao/projects_2026Q4/ WeQi_FTP/working_steps.sh，我这一轮还没读过。如果你想把这三条命令记进去，告诉我，我再动手。

sudo install -d -o gao -g sftp_users -m 755 /srv/external_data/client_sftp/files/Wei_Qi
sudo mount --bind /Work_bio/dropbox/Dropbox_Data/Wei_Qi /srv/external_data/client_sftp/files/Wei_Qi
sudo mount -o remount,bind,ro /srv/external_data/client_sftp/files/Wei_Qi

# 第三行把挂载设成只读，客户只能下载，改不了 Dropbox 里的原始数据。

这是临时挂载，重启后失效。要开机自动挂载，可以往 /etc/fstab 加一行 /Work_bio/dropbox/Dropbox_Data/Wei_Qi /srv/external_data/client_sftp/files/Wei_Qi none bind,ro 0 0，不过这属于系统改动，要你决定。

目前没有任何东西被复制，tmux 里也没有 rsync 在跑。
你执行完后告诉我，我会用 ls 和 md5sum 对照 MD5.csv 检查客户视角下的文件。

更正我上一条的一处说法：Quote_260428002 的挂载源是 .../Quote_260428002_mouse_rRNA//deleted，源目录已被删除，所以里面是空的（ls 显示 total 0）。我之前说"Wei Qi 能读到另一个客户的数据"过于严重。实际只剩一个空的残留挂载点，他最多看到这个目录名。不过 client_sftp 仍是所有客户共用的账号，这一点没变。

草稿箱里的邮件不用改，因为路径 /files/Wei_Qi 没变。


===
# working_steps.sh
# 有无关文件，请强制删除 /srv/external_data/client_sftp/files/Quote_260428002
# Bash Inspect stale mount before removal
# IN
# D=/srv/external_data/client_sftp/files/Quote_260428002; findmnt -o TARGET,SOURCE $D; ls -la $D; find $D -mindepth 1 | wc -l; ls -ld /home/gao/Dropbox/Quote_260428002* 2>&1 | head -3; grep -n Quote_260428002 /etc/fstab

# OUT
# TARGET                                               SOURCE
# /srv/external_data/client_sftp/files/Quote_260428002 /dev/sdb1[/dropbox/Dropbox_Data/Quote_260428002_mouse_rRNA//deleted]
# total 0
# bfs: error: /srv/external_data/client_sftp/files/Quote_260428002: No such file or directory.
# 0
# ls: cannot access '/home/gao/Dropbox/Quote_260428002*': No such file or directory
# 9:43 PM
# 这个目录是一个失效的 bind mount，不是数据，所以"删除"只需要卸载挂载点，再删掉空目录。

# 挂载源已不存在：源路径是 .../Quote_260428002_mouse_rRNA//deleted，find 报 No such file，里面 0 个文件。
# /home/gao/Dropbox 下也没有同名目录，不会误删任何数据。
# /etc/fstab 里没有这一项，卸载后重启不会再被挂回来。
# 需要你执行的 sudo 命令（我没有 sudo 权限）：


sudo umount -l /srv/external_data/client_sftp/files/Quote_260428002
sudo rmdir /srv/external_data/client_sftp/files/Quote_260428002

# 第一行用 -l（lazy）卸载，因为这个挂载源已经失效，普通卸载可能会报 busy。
# 第二行用 rmdir，而不是 rm -rf。它只能删空目录，万一里面还有东西会直接失败而不是删掉，所以比强制递归删除安全，也够用。
# 影响范围只是 files/ 下少一个空目录，不影响 Wei_Qi 的挂载，也不影响正在运行的 bore 隧道。
# 删完无法恢复，但它本来就是空的，没有可丢失的内容。
# 执行后告诉我，我再用 ls /srv/external_data/client_sftp/files 确认只剩 Wei_Qi。