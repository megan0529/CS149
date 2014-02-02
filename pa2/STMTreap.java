import java.util.concurrent.atomic.AtomicLong;

import org.deuce.Atomic;

public class STMTreap implements IntSet {
   static class Node {
      final int key;
      final int priority;
      Node      left;
      Node      right;
      
      Node(final int key, final int priority) {
         this.key = key;
         this.priority = priority;
      }
      
      public String toString() {
         return "Node[key=" + key + ", prio=" + priority + ", left=" + (left == null ? "null" : String.valueOf(left.key)) + ", right=" + (right == null ? "null" : String.valueOf(right.key)) + "]";
      }
   }
   
   private AtomicLong randState = new AtomicLong(0);
   private Node       root;
   
   @Atomic
   public boolean contains(final int key) {
      Node node = root;
      int nodeKey;
      while (node != null) {
         nodeKey = node.key;
         if (key == nodeKey) {
            return true;
         }
         node = key < nodeKey ? node.left : node.right;
      }
      return false;
   }
   
   @Atomic
   public void add(final int key) {
      root = addImpl(root, key);
   }
   
   @Atomic
   public void remove(final int key) {
      root = removeImpl(root, key);
   }
   
   @Atomic
   private Node addImpl(final Node node, final int key) {
      if (node == null) {
         return new Node(key, randPriority());
      }
      else {
         int nodeKey = node.key;
         if (key == nodeKey) {
            // no insert needed
            return node;
         }
         else if (key < nodeKey) {
            final Node leftNode = addImpl(node.left, key);
            if (node.left != leftNode) {
               node.left = leftNode;
               if (leftNode.priority > node.priority) {
                  return rotateRight(node);
               }
            }
            return node;
         }
         else {
            final Node rightNode = addImpl(node.right, key);
            if (node.right != rightNode) {
               node.right = rightNode;
               if (rightNode.priority > node.priority) {
                  return rotateLeft(node);
               }
            }
            return node;
         }
         
      }
   }
   
   private int randPriority() {
      // The constants in this 64-bit linear congruential random number
      // generator are from http://nuclear.llnl.gov/CNP/rng/rngman/node4.html
      while (true) {
         long old_rand = randState.get();
         long new_rand = old_rand * 2862933555777941757L + 3037000493L;
         if (randState.compareAndSet(old_rand, new_rand)) {
            return (int) (new_rand >> 30);
         }
      }
      
      // randState = randState * 2862933555777941757L + 3037000493L;
      // return (int) (randState >> 30);
   }
   
   @Atomic
   private Node rotateRight(final Node node) {
      // node nL
      // / \ / \
      // nL z ==> x node
      // / \ / \
      // x nLR nLR z
      final Node nL = node.left;
      node.left = nL.right;
      nL.right = node;
      return nL;
   }
   
   @Atomic
   private Node rotateLeft(final Node node) {
      final Node nR = node.right;
      node.right = nR.left;
      nR.left = node;
      return nR;
   }
   
   @Atomic
   private Node removeImpl(final Node node, final int key) {
      if (node == null) {
         // not present, nothing to do
         return null;
      }
      else {
         final Node oriLeftNode = node.left;
         final Node oriRightNode = node.right;
         final int nodeKey = node.key;
         
         if (key == nodeKey) {
            if (oriLeftNode == null) {
               return oriRightNode;
            }
            else if (oriRightNode == null) {
               return oriLeftNode;
            }
            else {
               // Two children, this is the hardest case. We will pretend
               // that node has -infinite priority, move it down, then retry
               // the removal.
               if (oriLeftNode.priority > oriRightNode.priority) {
                  // oriLeftNode needs to end up on top
                  final Node top = rotateRight(node);
                  final Node newTopRightNode = removeImpl(top.right, key);
                  if (top.right != newTopRightNode)
                     top.right = newTopRightNode;
                  return top;
               }
               else {
                  final Node top = rotateLeft(node);
                  final Node newTopLeftNode = removeImpl(top.left, key);
                  if (top.left != newTopLeftNode)
                     top.left = newTopLeftNode;
                  return top;
               }
            }
         }
         else if (key < nodeKey) {
            final Node leftNode = removeImpl(oriLeftNode, key);
            if (oriLeftNode != leftNode) {
               node.left = leftNode;
            }
            return node;
         }
         else {
            final Node rightNode = removeImpl(oriRightNode, key);
            if (oriRightNode != rightNode) {
               node.right = rightNode;
            }
            return node;
         }
      }
   }
}
