Artem Vasilyev  tema8
Song Han  songhan


We first created a ThreadPool class to schedule the 
accepted sockets to 8 threads. It's a waste to create 
a thread for each connection, but we are using these 8
threads to keep handling the connections.

The attributes of TheadPool are:
	private int                  numThreads;
    private LinkedList<Runnable> taskQueue;
    private WorkerThread[]       threads;

The prototype of ThreadPool's methods are:
	public ThreadPool(int _numThreads)  
	//constructor, spawning 8 threads and let them run
	
	public void submit(Runnable task)   
	//add a new taskRunnable to the task queue. Wake up
	//a thread that is waiting on it.

There's an innter class:
	private class WorkerThread extends Thread
	//waiting on the taskQueue to be non-empty, once 
	//waken up, pop the first task and let it run	


Here in order to let the task object access the handle()
method of ChatServer, the TaskRunnable class is made to
be an inner class of ChatServer. 

Synchronization:
1. The taskQueue is synchronized
2. The stateByName hashtable is synchronized
3. The history LinkedList is synchronized. Wait and Notify 
mechanism is also added to the addMessage() and 
recentMessages() methods


Bonus Point： (We implemented the 10 point version)

We assume that the first connection to 'all' room has to see
the history which meens that 'all room will be
created even if nobody is listneing to it. Overwise the first
user of 'all' room will not see messages prior
to his connection.


Our code archives 2nd level of Extra credit, because:
a) We use different locks for worker queue and ChatState, so
that requests are concurent and independent
b) Global lock is only used during post to 'all' room to 
prevent new rooms being created
c) Any post to normal room will only lock history of that 
room, all other rooms will be unlocked