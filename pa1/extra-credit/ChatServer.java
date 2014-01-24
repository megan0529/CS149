// ChatServer

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.net.ServerSocket;
import java.net.Socket;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChatServer {
   private static final Charset         utf8         = Charset.forName("UTF-8");
   
   private static final String          OK           = "200 OK";
   private static final String          NOT_FOUND    = "404 NOT FOUND";
   private static final String          HTML         = "text/html";
   private static final String          TEXT         = "text/plain";
   
   private static final Pattern         PAGE_REQUEST = Pattern.compile("GET /([^ /]+)/chat\\.html HTTP.*");
   private static final Pattern         PULL_REQUEST = Pattern.compile("POST /([^ /]+)/pull\\?last=([0-9]+) HTTP.*");
   private static final Pattern         PUSH_REQUEST = Pattern.compile("POST /([^ /]+)/push\\?msg=([^ ]*) HTTP.*");
   
   private static final String          CHAT_HTML;
   private static final int             NUM_THREAD   = 8;
   static {
      try {
         CHAT_HTML = getResourceAsString("chat.html");
      } catch (final IOException xx) {
         throw new Error("unable to start server", xx);
      }
   }
   
   private final int                    port;
   private final Map<String, ChatState> stateByName  = new HashMap<String, ChatState>();
   
   private ThreadPool                   threadPool;
   
   private class TaskRunnable implements Runnable {
      private Socket connSocket;
      
      public TaskRunnable(Socket _connSocket) {
         connSocket = _connSocket;
      }
      
      @Override
      public void run() {
         handle(connSocket);
      }
   };
   
   private class ThreadPool {
      private int                  numThreads;
      private LinkedList<Runnable> taskQueue;
      private WorkerThread[]       threads;
      
      public ThreadPool(int _numThreads) {
         numThreads = _numThreads;
         taskQueue = new LinkedList<Runnable>();
         threads = new WorkerThread[numThreads];
         for (int i = 0; i < numThreads; i++) {
            threads[i] = new WorkerThread();
            threads[i].start();
         }
      }
      
      public void submit(Runnable task) {
         synchronized (taskQueue) {
            taskQueue.addLast(task);
            taskQueue.notifyAll();
         }
      }
      
      // In order to let WorkerThread see the taskQueue, make the class inner
      // class
      private class WorkerThread extends Thread {
         @Override
         public void run() {
            // an infinite while loop keeping scrutinizing if the task queue is
            // empty. if not, pop out one task and let it run
            Runnable r;
            while (true) {
               synchronized (taskQueue) {
                  while (taskQueue.isEmpty()) {
                     try {
                        taskQueue.wait();
                     } catch (InterruptedException e) {
                        e.printStackTrace();
                     }
                  }
                  r = taskQueue.pollFirst();
               }
               r.run();
            }
            
         }
      };
   };
   
   /**
    * Constructs a new {@link ChatServer} that will service requests on the
    * specified <code>port</code>. <code>state</code> will be used to hold the
    * current state of the chat.
    */
   public ChatServer(final int port) throws IOException {
      this.port = port;
      threadPool = new ThreadPool(NUM_THREAD);
   }
   
   public void runForever() throws Exception {
      final ServerSocket server = new ServerSocket(port);
      while (true) {
         final Socket connSocket = server.accept();
         TaskRunnable t = new TaskRunnable(connSocket);
         threadPool.submit(t);
      }
   }
   
   private void handle(final Socket connection) {
      try {
         final BufferedReader xi = new BufferedReader(new InputStreamReader(connection.getInputStream()));
         final OutputStream xo = connection.getOutputStream();
         
         final String request = xi.readLine();
         System.out.println(Thread.currentThread() + ": " + request);
         
         Matcher m;
         if (request == null) {
            sendResponse(xo, NOT_FOUND, TEXT, "Empty request.");
         }
         else if (PAGE_REQUEST.matcher(request).matches()) {
            sendResponse(xo, OK, HTML, CHAT_HTML);
         }
         
         // PULL
         else if ((m = PULL_REQUEST.matcher(request)).matches()) {
            final String room = m.group(1);
            final long last = Long.valueOf(m.group(2));
            sendResponse(xo, OK, TEXT, getState(room).recentMessages(last));
         }
         // PUSH
         else if ((m = PUSH_REQUEST.matcher(request)).matches()) {
            final String room = m.group(1);
            final String msg = m.group(2);
            
            // if (room.equals("all")) {
            // getState(room).addMessage(msg);
            // synchronized (stateByName) {
            // for (String key: stateByName.keySet()) {
            // if(key != "all") {
            // getState(key).addMessage(msg);
            // }
            // }
            // }
            // else {
            // getState(room).addMessage(msg);
            // getState("all").addMessage(msg);
            // }
            
            getState(room).addMessage(msg);
            sendResponse(xo, OK, TEXT, "ack");
         }
         else {
            sendResponse(xo, NOT_FOUND, TEXT, "Malformed request.");
         }
         
         connection.close();
      } catch (final Exception xx) {
         xx.printStackTrace();
         try {
            connection.close();
         } catch (final Exception yy) {
            // Silently discard exceptions while closing the socket.
         }
      }
   }
   
   /**
    * Writes a minimal but valid HTML response to <code>output</code>.
    */
   private static void sendResponse(final OutputStream output, final String status, final String contentType, final String content) throws IOException {
      final byte[] data = content.getBytes(utf8);
      final String headers = "HTTP/1.0 " + status + "\n" + "Content-Type: " + contentType + "; charset=utf-8\n" + "Content-Length: " + data.length + "\n\n";
      
      final BufferedOutputStream xo = new BufferedOutputStream(output);
      xo.write(headers.getBytes(utf8));
      xo.write(data);
      xo.flush();
      
      System.out.println(Thread.currentThread() + ": replied with " + data.length + " bytes");
   }
   
   private ChatState getState(final String room) {
      System.out.println("room is " + room);
      synchronized (stateByName) {
         ChatState state = stateByName.get(room);
         if (state == null) {
            state = new ChatState(room);
            stateByName.put(room, state);
         }
         return state;
      }
   }
   
   /**
    * Reads the resource with the specified name as a string, and then returns
    * the string. Resource files are searched for using the same classpath
    * mechanism used for <code>.class</code> files, so they can either be in the
    * same directory as bare <code>.class</code> files or included in the
    * <code>.jar</code> file.
    */
   private static String getResourceAsString(final String name) throws IOException {
      final Reader xi = new InputStreamReader(ChatServer.class.getClassLoader().getResourceAsStream(name));
      try {
         final StringBuffer result = new StringBuffer();
         final char[] buf = new char[8192];
         int n;
         while ((n = xi.read(buf)) > 0) {
            result.append(buf, 0, n);
         }
         return result.toString();
      } finally {
         try {
            xi.close();
         } catch (final IOException xx) {
            // Discard exceptions.
         }
      }
   }
   
   /**
    * Runs a chat server, with a default port of 8080.
    */
   public static void main(final String[] args) throws Exception {
      final int port = args.length == 0 ? 8080 : Integer.parseInt(args[0]);
      new ChatServer(port).runForever();
   }
}
