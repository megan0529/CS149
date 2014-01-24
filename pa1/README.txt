Artem Vasilyev  tema8
Song Han


We assume that the first connection to 'all' room has to see the history which meens that 'all room will be
created even if nobody is listneing to it. Overwise the first user of 'all' room will not see messages prior
to his connection.


Our code archives 2nd level of Extra credit, because:
a) We use different locks for workes que and ChatState, so that requests are concurent and independent
b) Global lock is only used during post to 'all' room to prevent new rooms being created
c) Any post to normal room will only lock history of that room, all other rooms will be unlocked