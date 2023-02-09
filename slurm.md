## SLURM for Debian

### Munge

SLURM is the name of the scheduling software used to queue and manage the jobs we submit to our clusters.
You can either build it from source or get the appropriate packages from the package manager of you distribution
(apt for debian, yum for rocky)
As a word of warning: most of these commands will require root (sudo) privileges.
The basic guideline is detailed in e.g. https://gist.github.com/asmateus/301b0cb86700cbe74c269b27f2ecfbef

First, install `munge` (this is used for authentication) on the cluster head:
```
sudo apt-get install libmunge-dev libmunge2 munge
```

At this point a `munge` user should have been created.
Then generate a key that `munge` will use to authenticate:
```
sudo dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key
```
Set the permission of the appropriate directory so that munge has access to its key:
```
sudo chown -R /etc/munge/
sudo chmod 400 /etc/munge/munge.key
```

This key should be copied to every node on the cluster:
```
scp /etc/munge/munge.key <nodename>:/etc/munge/munge.key
sudo ssh <nodename> 'chown -R /etc/munge/ ; chmod 400 /etc/munge/munge.key'
```

The `munge` account should not have login access, so the passwd file has to be edited (
and all nodes):
```
sudo vim /etc/passwd
munge:x:<num>:<num>::/nonexistent:/bin/false
```

Then the `munge` daemon can be started on every node, including the head:
```
sudo systemctl enable munge
sudo systemctl start munge
```

If everything goes smoothly, the following commands should return a success:
```
munge -n | unmunge
munge -n | ssh <nodename> unmunge
```

If not, you should look into the logs (either via `systemctl status munge` or at /var/log/munge/munged.log).
Typical reasons for the service not starting may include missing permissions for the `munge` user for certain directories.
Should that be the case, change the ownership of said directories or files via `sudo chown muge:munge /path/to/directory`

### Slurm

Next, the `slurm` package should be installed:
```
sudo apt-get install slurm-wlm
```

This should create a `slurm` user as well as configuration (`/etc/slurm-llnl`) and log (`/var/log/slurm-llnl`) directories.
If no `slurm.conf` is created under `/etc/slurm-llnl`, you have to create a new file. Use either templates from online sources
or one from our existing older cluster heads. The configuration and log directories should be owned by the `slurm` user, as well
as any other files that `slurm` needs to access.

In order to create a PID for `slurm`, you have to create the files on the head
```
sudo touch /var/run/slurmctld.pid
sudo chown slurm:slurm /var/run/slurmctld.pid
```
or whichever directory you have specified in the `slurm.conf` file.

On the compute nodes, you have to create a pid for `slurmd`
```
sudo touch /var/run/slurmd.pid
sudo chown slurm:slurm /var/run/slurmd.pid
```

Generally, the `slurmctld` daemon should only ever interact with the cluster head, and the `slurmd` daemon with the compute nodes.
If you want to use the head not just for submission but also for computing jobs, you could also have it running on the head,
but the gain is usually outweighed by the complications that creates.

*Note: the linked tutorial tells to specify the `slurm` user in the .services files (`/lib/systemd/system/slurmd.service` or
`slurmctld.service`), assigning different users to the nodes and the head. It is also possible to specify the `SlurmUser` in the
`slurm.conf` file, and just have slurm handle everything on the head and the nodes.*

If you plan on using not just CPUs but GPUs as well, you will have to tell slurm to look for these generic resources.
On each compute node, create a file under `/etc/slurm-llnl/`
```
sudo vim /etc/slurm-llnl/gres.conf
# enter the generic resource specifications. For NVIDIA GPUs, that would be
Name=gpu File=/dev/nvidia0
# For multiple GPU systems, you can specify each graphics card as
# File=/dev/nvidia1 etc.
```

Once that is taken care of, you can start the services. On the head that includes
```
sudo systemctl enable slurmctld
sudo systemctl start slurmctld
```

And for the nodes
```
sudo systemctl enable slurmd
sudo systemctl start slurmd
```

You should check with the `systemctl status` command whether the services are up and running. If not, look at the log files specified
in `slurm.conf`, or try running the daemon directly as the `slurm` user:
```
sudo -u slurm slurmctld -Dvvv
```
This will generate some logs that are generally easier to interpret than the `systemctl` messages. Increasing the number of `v`s increases
the verbosity of the output.

If everyhting is up and running, you can try submitting simple jobs to a node to see whether the cluster is configured properly.

### SlurmDBD

At this point, you should have a working slurm installation. If you wish to take advantage of the advanced scheduling algorithm
(enabled by `PriorityType=priority/multifactor` in `slurm.conf`), you will have to create a database and install a third daemon
on the cluster head.
This guide largely follows https://github.com/nateGeorge/slurm_gpu_ubuntu, with some caveats.

Slurm needs a database to compute some parameters, such as job age in order to fully leverage its capabilities.
If they are not already installed, then you should get `mysql` and `mariadb`
```
sudo apt-get install mariadb-client mariadb-server
```
Enable and start the services:
```
sudo systemctl enable mariadb
sudo systemctl start mariadb
```

If everything is running smoothly, the database can be created for slurm:
```
sudo mysql -u root
# Create the database and add the slurm user ; the > character is the mysql prompt
> create database slurm_acct_db;
> create user 'slurm'@'localhost';
# If you want to set a password, you will have to include the password in slurmdbd.conf
# In that case, slurmdbd.conf should not be readable to anyone other than the slurm user
# Alternatively, you could also skip the password and live with the consequences
set password for 'slurm'@'localhost' = password('slurmdbpass');
# give slurm access to the things it requires
> grant usage on *.* to 'slurm'@'localhost';
> grant all privileges on slurm_acct_db.* to 'slurm'@'localhost';
> flush privileges;
> exit
```

Now you will need to create `slurmdbd.conf` in `/etc/slurm-llnl/`. You can either use the template of
the linked tutorial, or one from an older cluster head. The same thing also applies for `slurmdbd.services`,
which should be placed under `/lib/systemd/system/`.

Make sure that the files and directories referred to in the configuration and services exist (such as the log
and pid files) and that they are owned by the `slurm` user. If you have set a password for the the `slurmdbd`
during the mysql session, make sure you include that in `slurmdbd.conf` as `StoragePass=slurmdbdpass`, or whatever
password you chose.

At this point you can start the `slurmdbd` daemon
```
sudo systemctl enable slurmdbd
sudo systemctl start slurmdbd
```
Check with `sudo systemctl status slurmdbd` whether it is actually running (sometimes `systemctl` does not
give any error upon start, even though it fails). The same general steps of troubleshooting apply as with
the other daemons: look in the log files (`/var/log/slurm-llnl/slurmdbd.log`), or start the service
directly as the `slurm` user with `sudo -u slurm slurmdbd -Dvvv`.

Assuming everything is up and running, you can restart the slurm services. On every compute node, do
```
sudo systemctl restart slurmd
```
On the head
```
sudo systemctl restart slurmctld
sudo systemctl restart slurmdbd
```

### Configuring the cluster

Now it is time to add the cluster to the slurm database. An exhaustive guide can be found at e.g.
https://wiki.fysik.dtu.dk/niflheim/Slurm_accounting

The first step is to define our cluster:
```
sacctmgr add cluster clustername
```
`clustername` should  be the same as in the `slurm.conf` and `/var/lib/slurm-llnl/slurmctld/clustername`

Then add an account (a group that users belong to)
```
sacctmgr add account mdyusers Description="mdyusers" Organization=mdyusers
```

Finally, add everyone to the account and cluster
```
sacctmgr create user name=<username> cluster=clustername DefaultAccount=mdyusers Account=mdyusers
```

Once everything in in place, the queueing should now be enabled. You may need to resatrt the daemons
(`slurmd` on the compute nodes as well as `slurmctld` and `slurmdbd` on the head) before the changes take effect.

### Troubleshooting remarks

If nodes are put in a `down` state without any apparent reason, but they can be reached via ping and ssh,
check whether the system times are in sync with the clusterhead (that is, compare the output of
```
date
```
on the head and the node in question). If those are out of sync by more than 5 minutes, munge won't allow any
communication between the controller and slurmd. Set the hwclock time manually by entering
```
date --set hh:mm:ss
```
where hh:mm:ss is the current time as reported by the head.
